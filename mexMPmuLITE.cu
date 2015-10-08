/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. B=mexFunction(A).
 *
 * Copyright 2012 The MathWorks, Inc.
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cstdlib>
#include <algorithm>
#include <iostream>
using namespace std;

const int nt0 = 61,  Nthreads = 1024, lockout = nt0-1, NchanMax = 128, block = 32, NrankMax = 3;
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	Conv1D(const double *Params, const float *data, const float *W, float *conv_sig){    
  volatile __shared__ float  sW[nt0*NrankMax], sdata[(Nthreads+nt0)*NrankMax]; 
  float x;
  int tid, tid0, bid, i, nid, Nrank, NT, Nfilt;

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  Nfilt    	=   (int) Params[1];
  NT      	=   (int) Params[0];
  Nrank     = (int) Params[6];
  
  if(tid<nt0*((int) Params[6]))
      sW[tid]= W[tid%nt0 + (bid + Nfilt * (tid/nt0))* nt0];
  __syncthreads();
  
  tid0 = 0;
  while (tid0<NT-Nthreads-nt0+1){
	  if (tid<nt0*NrankMax) sdata[tid%nt0 + (tid/nt0)*(Nthreads+nt0)] = 
			data[tid0 + tid%nt0+ NT*(bid + Nfilt*(tid/nt0))];
	  #pragma unroll 3
      for(nid=0;nid<Nrank;nid++){
          sdata[tid + nt0+nid*(Nthreads+nt0)] = data[nt0+tid0 + tid+ NT*(bid +nid*Nfilt)];
	  }
	  __syncthreads();
      
	  x = 0.0f;
      for(nid=0;nid<Nrank;nid++){
		  #pragma unroll 4
          for(i=0;i<nt0;i++)
              x    += sW[i + nid*nt0] * sdata[i+tid + nid*(Nthreads+nt0)];
	  }
      conv_sig[tid0  + tid + NT*bid]   = x;
      
      tid0+=Nthreads;
      __syncthreads();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void  bestFilter(const double *Params, const float *data, 
	const float *mu, const float *lam, float *xbest, float *err, int *ftype){
  int tid, tid0, i, bid, NT, Nfilt, ibest = 0;
  float Th,  Cf, Ci, xb, Cbest = 0.0f;

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  NT 		= (int) Params[0];
  Nfilt 	= (int) Params[1];
  Th 		= (float) Params[2];
 
  tid0 = tid + bid * Nthreads;
  if (tid0<NT){
    for (i=0; i<Nfilt;i++){
      Ci = data[tid0 + NT * i] + mu[i] * lam[i];
      Cf = Ci * Ci / (lam[i] + 1.0f) - lam[i]*mu[i]*mu[i];
		if (Cf > Cbest){
			Cbest 	= Cf;
			xb      = Ci - lam[i] * mu[i]; // /(lam[i] + 1);
			ibest 	= i;
		}
    }
    if (Cbest > Th*Th){
      err[tid0] 	= Cbest;
      xbest[tid0] 	= xb;
      ftype[tid0] 	= ibest;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	cleanup_spikes(const double *Params, const float *xbest, const float *err, 
	const int *ftype, int *st, int *id, float *x, float *C, int *counter){
  int indx, maxFR, NTOT, tid, bid, NT, tid0,  j;
  volatile __shared__ float sdata[Nthreads+2*lockout+1];
  bool flag=0;
  float err0;
  
  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  
  NT      	=   (int) Params[0];
  maxFR 	= (int) Params[3];
  tid0 		= bid * Nthreads;


  if(tid0<NT-Nthreads-2*lockout-1){       
    if (tid<2*lockout)
      sdata[tid] = err[tid0 + tid];
    sdata[tid+2*lockout] = err[2*lockout + tid0 + tid];

    __syncthreads();
    
    err0 = sdata[tid+lockout];
    if(err0>1e-10){
        flag = 0;
        for(j=-lockout;j<=lockout;j++)
            if(sdata[tid+lockout+j]>err0){
                flag = 1;
                break;
            }
        if(flag==0){
            indx = atomicAdd(&counter[0], 1);
            if (indx<maxFR){
                st[indx] = tid+lockout         + tid0;
                id[indx] = ftype[tid+lockout   + tid0];
                x[indx]  = xbest[tid+lockout     + tid0];
                C[indx]  = err0;
            }
        }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	subSpikes(const double *Params, const int *st, const int *id, const float *x, const int *counter, float *dout, const float *WtW){
  int tid, bid,  NT, ind, tcurr, Nfilt;
  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  NT 		= (int) Params[0];
  Nfilt 	= (int) Params[1];

  for(ind=counter[1]; ind<counter[0];ind++){
    tcurr = tid + st[ind]-nt0+1;
    if (tcurr>=0 & tcurr<NT)
      dout[tcurr + bid*NT] -= x[ind] * WtW[tid + id[ind]*(2*nt0-1) + (2*nt0-1)*Nfilt*bid];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	subtract_spikes(const double *Params,  const int *st, const int *id, const float *x, const int *counter, float *dataraw, const float *W, const float *U){
  int tid, bid, Nblocks, i, NT, ind, Nchan;
  __shared__ float sh_W[nt0], sh_U[NchanMax];
  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  Nblocks   = gridDim.x;
  NT        = (int) Params[0];
  Nchan     = (int) Params[5];
  ind       = bid;

  while(ind<counter[0]){
    if (tid<nt0) sh_W[tid] = W[tid + nt0*id[ind]];
    sh_U[tid] = U[tid + Nchan*id[ind]];

    __syncthreads();
    for (i=0;i<nt0;i++)
      dataraw[i + st[ind] + NT * tid] -= x[ind] * sh_W[i] * sh_U[tid];
    ind+= Nblocks;
    __syncthreads();
  }

}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void getWgradient(const double *Params, const int *st, const int *id, const float *x,  const int *counter, const float *datarez, const float *U, float *dW){
  int tid, bid, i, ind, NT, Nchan;
  float xprod; 
  volatile __shared__ float sh_U[NchanMax];
  NT = (int) Params[0];
    Nchan = (int) Params[5];

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  while(tid<Nchan){
    sh_U[tid] = U[tid + bid*Nchan];
    tid+= blockDim.x;
  }
  tid 		= threadIdx.x;
  __syncthreads();

  for(ind=0; ind<counter[0];ind++)
    if (id[ind]==bid){
      xprod = 0.0f;
      for (i=0;i<Nchan;i++)
	xprod+= sh_U[i] * datarez[st[ind] + tid + NT * i];
      dW[tid + nt0 * bid] += xprod * x[ind];
    }
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void getUgradient(const double *Params, const int *st, const int *id, const float *x,  const int *counter, const float *datarez, const float *W, float *dU){  
  int j, tid, bid, i, ind, NT, Nchan;
  float xprod; 
  volatile __shared__ float sh_M[NchanMax*nt0], sh_W[nt0];

  NT = (int) Params[0];
    Nchan = (int) Params[5];

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  if (tid<nt0)
    sh_W[tid] = W[tid + nt0*bid];
 
  __syncthreads();

  for(ind=0; ind<counter[0];ind++)
    if (id[ind]==bid){
      if(tid<nt0)
	for (j=0;j<Nchan;j++)
	  sh_M[tid + nt0*j] = datarez[tid + st[ind] + NT*j];
      __syncthreads();

      xprod = 0.0f;
      for (i=0;i<nt0;i++)
	xprod+= sh_W[i] * sh_M[i + tid*nt0];
      dU[tid + bid*Nchan] += xprod * x[ind];
      __syncthreads();
    }  
}
//////////////////////////////////////////////////////////////////////////////////////////

/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
  /* Declare input variables*/
  double *Params, *d_Params;
  int blocksPerGrid, NT, maxFR, Nchan;
  int const threadsPerBlock = Nthreads;

  /* Initialize the MathWorks GPU API. */
  mxInitGPU();

  /* read Params and copy to GPU */
  Params  	= (double*) mxGetData(prhs[0]);
  NT		= (int) Params[0];
  blocksPerGrid	= (int) Params[1];
  maxFR		= (int) Params[3];
   Nchan = (int) Params[5];

  cudaMalloc(&d_Params,      sizeof(double)*mxGetNumberOfElements(prhs[0]));
  cudaMemcpy(d_Params,Params,sizeof(double)*mxGetNumberOfElements(prhs[0]),cudaMemcpyHostToDevice);

  /* collect input GPU variables*/
  mxGPUArray const  *W,   *data,   *WtW, *mu,   *lam;
  const float     *d_W, *d_data, *d_WtW,  *d_mu, *d_lam;
  
  data       = mxGPUCreateFromMxArray(prhs[1]);
  d_data     = (float const *)(mxGPUGetDataReadOnly(data));
  W             = mxGPUCreateFromMxArray(prhs[2]);
  d_W        	= (float const *)(mxGPUGetDataReadOnly(W));
  WtW       	= mxGPUCreateFromMxArray(prhs[3]);
  d_WtW     	= (float const *)(mxGPUGetDataReadOnly(WtW));
  mu            = mxGPUCreateFromMxArray(prhs[4]);
  d_mu          = (float const *)(mxGPUGetDataReadOnly(mu));
  lam       	= mxGPUCreateFromMxArray(prhs[5]);
  d_lam     	= (float const *)(mxGPUGetDataReadOnly(lam));

  /* allocate new GPU variables*/  
  float *d_err,*d_C, *d_xbest, *d_x, *d_dout;
  int *d_st,  *d_ftype,  *d_id, *d_counter;

  cudaMalloc(&d_dout,   NT * blocksPerGrid* sizeof(float));

  cudaMalloc(&d_err,   NT * sizeof(float));
  cudaMalloc(&d_xbest,   NT * sizeof(float));
  cudaMalloc(&d_ftype, NT * sizeof(int));
  cudaMalloc(&d_st,    maxFR * sizeof(int));
  cudaMalloc(&d_id,    maxFR * sizeof(int));
  cudaMalloc(&d_x,     maxFR * sizeof(float));
  cudaMalloc(&d_C,     maxFR * sizeof(float));
  cudaMalloc(&d_counter,   2*sizeof(int));
 
  cudaMemset(d_dout,    0, NT * blocksPerGrid * sizeof(float));
  cudaMemset(d_counter, 0, 2*sizeof(int));
  cudaMemset(d_st,      0, maxFR *   sizeof(int));
  cudaMemset(d_id,      0, maxFR *   sizeof(int));
  cudaMemset(d_x,       0, maxFR *    sizeof(float));
  cudaMemset(d_C,       0, maxFR *    sizeof(float));

  int *counter;
  counter = (int*) calloc(1,sizeof(int));
 
  Conv1D<<<blocksPerGrid,threadsPerBlock>>>(d_Params, d_data, d_W, d_dout); 
  for(int k=0;k<(int) Params[4];k++){
    cudaMemset(d_err,     0, NT * sizeof(float));
    cudaMemset(d_ftype,   0, NT * sizeof(int));
    cudaMemset(d_xbest,   0, NT * sizeof(float));

    bestFilter<<<NT/Nthreads,threadsPerBlock>>>(    d_Params, d_dout, d_mu, d_lam, d_xbest, d_err, d_ftype);
    cleanup_spikes<<<NT/Nthreads,threadsPerBlock>>>(d_Params, d_xbest, d_err, d_ftype, d_st, d_id, d_x, d_C, d_counter);
 
    cudaMemcpy(counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    if (counter[0]>maxFR){
      counter[0] = maxFR;
      cudaMemcpy(d_counter, counter, sizeof(int), cudaMemcpyHostToDevice);      
    }
    
    subSpikes<<<blocksPerGrid, 2*nt0-1>>>(d_Params, d_st, d_id, d_x, d_counter, d_dout,    d_WtW);

    cudaMemcpy(d_counter+1, d_counter, sizeof(int), cudaMemcpyDeviceToHost);

    if(counter[0]==maxFR)
      break;
  }

  float *x, *C;
  int *st, *id;
  int minSize;
  if (counter[0]<maxFR)  minSize = counter[0];
  else                   minSize = maxFR;
  const mwSize dimst[] 	= {minSize,1}; 
  plhs[0] = mxCreateNumericArray(2, dimst, mxINT32_CLASS, mxREAL);
  st = (int*) mxGetData(plhs[0]);
  plhs[1] = mxCreateNumericArray(2, dimst, mxINT32_CLASS, mxREAL);
  id = (int*) mxGetData(plhs[1]);
  plhs[2] = mxCreateNumericArray(2, dimst, mxSINGLE_CLASS, mxREAL);
  x =  (float*) mxGetData(plhs[2]);
  plhs[3] = mxCreateNumericArray(2, dimst, mxSINGLE_CLASS, mxREAL);
  C =  (float*) mxGetData(plhs[3]);
  cudaMemcpy(st, d_st, minSize * sizeof(int),   cudaMemcpyDeviceToHost);
  cudaMemcpy(id, d_id, minSize * sizeof(int),   cudaMemcpyDeviceToHost);
  cudaMemcpy(x,   d_x, minSize * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(C,   d_C, minSize * sizeof(float), cudaMemcpyDeviceToHost);

  cudaFree(d_ftype);
  cudaFree(d_err);
  cudaFree(d_xbest);
  cudaFree(d_st);
  cudaFree(d_id);
  cudaFree(d_x);
  cudaFree(d_C);
  cudaFree(d_counter);
  cudaFree(d_Params);

  cudaFree(d_dout);

  mxGPUDestroyGPUArray(data);
  mxGPUDestroyGPUArray(WtW);
  mxGPUDestroyGPUArray(W);
  mxGPUDestroyGPUArray(mu);
  mxGPUDestroyGPUArray(lam);
}
