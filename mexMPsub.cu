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

const int nt0 = 61,  Nthreads = 1024,   lockout = nt0-1, NchanMax = 128, block = 32;

//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	Conv1D(const double *Params, const float *data, const float *W, float *conv_sig){    
  __shared__ float sW[nt0], sdata[Nthreads+nt0]; 
  float x;
  int tid, tid0, bid, i, NT;

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  
  if(tid<nt0)        sW[tid]= W[tid + bid * nt0];
  __syncthreads();
	 	 
  NT      	=   (int) Params[0];
  tid0 = 0;
  while (tid0<NT-Nthreads-nt0+1){    
    if (tid<nt0) sdata[tid] = data[tid0 + tid+ NT*bid];
    sdata[nt0+tid] = data[nt0+tid0 + tid+ NT*bid];
     __syncthreads();
    
    x = 0.0f;
    for(i=0;i<nt0;i++)
      x    += sW[i] * sdata[i+tid];
    
    conv_sig[tid0  + tid + NT*bid]   = x;
    
    tid0+=Nthreads;
     __syncthreads();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void  bestFilter(const double *Params, const float *data, float *err, int *ftype){

  int tid, tid0, i, bid, NT, Nfilt, ibest = 0;
  float xbest = 0.0f, Th;

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  NT 		= (int) Params[0];
  Nfilt 	= (int) Params[1];
  Th 		= (float) Params[2];

  tid0 = tid + bid * Nthreads;
  if (tid0<NT){
    for (i=0; i<Nfilt;i++)
      if (abs(data[tid0 + NT * i]) > abs(xbest)){
	xbest = data[tid0 + NT * i];
	ibest = i;
      }
    if (abs(xbest)>Th){
      err[tid0] 	= xbest;
      ftype[tid0] 	= ibest;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	cleanup_spikes(const double *Params, const float *err, const int *ftype, int *st, int *id, float *x, int *counter){
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
      sdata[tid] = abs(err[tid0 + tid]*err[tid0 + tid]);
    sdata[tid+2*lockout] = abs(err[2*lockout + tid0 + tid]*err[2*lockout + tid0 + tid]);

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
	    x[indx]  = err[tid+lockout     + tid0];
	  }
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void	subSpikes(const double *Params, const int *st, const int *id, const float *x, const int *counter, float *dout, const float *WtW){
  int tid, bid,  NT, ind, tcurr, Nfilt, Nchan;
  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  NT 		= (int) Params[0];
  Nfilt 	= (int) Params[1];
  Nchan         = (int) Params[5];

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
  Nblocks       = gridDim.x;
  NT = (int) Params[0];
  Nchan         = (int) Params[5];
  ind = bid + counter[1];

  while(ind<counter[0]){
    while (tid<nt0){ sh_W[tid] = W[tid + nt0*id[ind]]; tid+=blockDim.x;}    
    tid 		= threadIdx.x;
    sh_U[tid] = U[tid + Nchan*id[ind]];
    
    __syncthreads();
    for (i=0;i<nt0;i++)
      dataraw[i + st[ind] + NT * tid] -= x[ind] * sh_W[i] * sh_U[tid];
    ind+= Nblocks;
    __syncthreads();
  }

}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void getWgradient(const double *Params, const int *st, const int *id, 
        const float *x,  const int *counter, const float *datarez, const float *U, float *dW){
  int tid, bid, i, ind, NT, Nchan;
  float xprod; 
  volatile __shared__ float sh_U[NchanMax];
  NT = (int) Params[0];
    Nchan         = (int) Params[5];

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
  Nchan         = (int) Params[5];

  tid 		= threadIdx.x;
  bid 		= blockIdx.x;
  while(tid<nt0){
    sh_W[tid] = W[tid + nt0*bid];
    tid+=blockDim.x;
  }
  tid 		= threadIdx.x;
 
  __syncthreads();
  
  for(ind=0; ind<counter[0];ind++)
      if (id[ind]==bid){
          while(tid<nt0){
              for (j=0;j<Nchan;j++)
                  sh_M[tid + nt0*j] = datarez[tid + st[ind] + NT*j];
              tid+=blockDim.x;
          }
          tid 		= threadIdx.x;
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
  Params        = (double*) mxGetData(prhs[0]);
  NT            = (int) Params[0];
  blocksPerGrid	= (int) Params[1];
  maxFR         = (int) Params[3];
  Nchan         = (int) Params[5];
  cudaMalloc(&d_Params,      sizeof(double)*mxGetNumberOfElements(prhs[0]));
  cudaMemcpy(d_Params,Params,sizeof(double)*mxGetNumberOfElements(prhs[0]),cudaMemcpyHostToDevice);
  
  /* collect input GPU variables*/
  mxGPUArray const  *W,   *U,   *dataraw,   *data,   *WtW;
  const float     *d_W, *d_U, *d_dataraw, *d_data, *d_WtW;
  
  dataraw       = mxGPUCreateFromMxArray(prhs[1]);
  d_dataraw     = (float const *)(mxGPUGetDataReadOnly(dataraw));
  W             = mxGPUCreateFromMxArray(prhs[2]);
  d_W        	= (float const *)(mxGPUGetDataReadOnly(W));
  U         	= mxGPUCreateFromMxArray(prhs[3]);
  d_U        	= (float const *)(mxGPUGetDataReadOnly(U));
  data        	= mxGPUCreateFromMxArray(prhs[4]);
  d_data        = (float const *)(mxGPUGetDataReadOnly(data));
  WtW       	= mxGPUCreateFromMxArray(prhs[5]);
  d_WtW     	= (float const *)(mxGPUGetDataReadOnly(WtW));
  
  /* allocate new GPU variables*/
  float *d_err, *d_x, *d_dout;
  int *d_st, *d_ftype,  *d_id, *d_counter;
  
  cudaMalloc(&d_dout,   NT * blocksPerGrid* sizeof(float));

  cudaMalloc(&d_err,   NT * sizeof(float));
  cudaMalloc(&d_ftype, NT * sizeof(int));
  cudaMalloc(&d_st,    maxFR * sizeof(int));
  cudaMalloc(&d_id,    maxFR * sizeof(int));
  cudaMalloc(&d_x,     maxFR * sizeof(float));
  cudaMalloc(&d_counter,   2*sizeof(int));
 
  cudaMemset(d_dout,    0, NT * blocksPerGrid * sizeof(float));
  cudaMemset(d_counter, 0, 2*sizeof(int));
  cudaMemset(d_st,      0, maxFR *   sizeof(int));
  cudaMemset(d_id,      0, maxFR *   sizeof(int));
  cudaMemset(d_x,       0, maxFR *    sizeof(float));


  mxGPUArray *datarez, *dW, *dU;
  float *d_datarez, *d_dW, *d_dU;
  const mwSize dimsu[] 	= {NT,Nchan}; 
  datarez 		= mxGPUCreateGPUArray(2, dimsu, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);  
  d_datarez 		= (float *)(mxGPUGetData(datarez));
  cudaMemcpy(d_datarez, d_dataraw,  NT * Nchan * sizeof(float), cudaMemcpyDeviceToDevice);

  const mwSize dimsdW[] = {nt0,blocksPerGrid}; 
  dW 		= mxGPUCreateGPUArray(2, dimsdW, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);  
  d_dW 		= (float *)(mxGPUGetData(dW));
  cudaMemset(d_dW, 0,  nt0*blocksPerGrid * sizeof(float));

  const mwSize dimsdU[] = {Nchan,blocksPerGrid}; 
  dU 		= mxGPUCreateGPUArray(2, dimsdU, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);  
  d_dU 		= (float *)(mxGPUGetData(dU));
  cudaMemset(d_dU, 0,  Nchan*blocksPerGrid * sizeof(float));

  int *counter;
  counter = (int*) calloc(1,sizeof(int));
 
  Conv1D<<<blocksPerGrid,threadsPerBlock>>>(d_Params, d_data, d_W, d_dout); 
  for(int k=0;k<(int) Params[4];k++){
    cudaMemset(d_err,     0, NT * sizeof(float));
    cudaMemset(d_ftype,   0, NT * sizeof(int));

    bestFilter<<<NT/Nthreads,threadsPerBlock>>>(    d_Params, d_dout, d_err, d_ftype);
    cleanup_spikes<<<NT/Nthreads,threadsPerBlock>>>(d_Params, d_err, d_ftype, d_st, d_id, d_x, d_counter);
 
    cudaMemcpy(counter, d_counter, sizeof(int), cudaMemcpyDeviceToHost);
    if (counter[0]>maxFR){
      counter[0] = maxFR;
      cudaMemcpy(d_counter, counter, sizeof(int), cudaMemcpyHostToDevice);      
    }
      
    subtract_spikes<<<128,Nchan>>>(       d_Params, d_st, d_id, d_x, d_counter, d_datarez, d_W, d_U);
    subSpikes<<<blocksPerGrid, 2*nt0-1>>>(d_Params, d_st, d_id, d_x, d_counter, d_dout,    d_WtW);

    cudaMemcpy(d_counter+1, d_counter, sizeof(int), cudaMemcpyDeviceToHost);

    if(counter[0]==maxFR)
      break;
  }

  getWgradient<<<blocksPerGrid,nt0>>>(  d_Params, d_st, d_id, d_x, d_counter, d_datarez, d_U, d_dW);
  getUgradient<<<blocksPerGrid,Nchan>>>(d_Params, d_st, d_id, d_x, d_counter, d_datarez, d_W, d_dU);

  plhs[0] 	= mxGPUCreateMxArrayOnGPU(datarez);
  plhs[1] 	= mxGPUCreateMxArrayOnGPU(dW);
  plhs[2] 	= mxGPUCreateMxArrayOnGPU(dU);


  float *x;
  int *st, *id;
  int minSize;
  if (counter[0]<maxFR)  minSize = counter[0];
  else                   minSize = maxFR;
  const mwSize dimst[] 	= {minSize,1}; 
  plhs[3] = mxCreateNumericArray(2, dimst, mxINT32_CLASS, mxREAL);
  st = (int*) mxGetData(plhs[3]);
  plhs[4] = mxCreateNumericArray(2, dimst, mxINT32_CLASS, mxREAL);
  id = (int*) mxGetData(plhs[4]);
  plhs[5] = mxCreateNumericArray(2, dimst, mxSINGLE_CLASS, mxREAL);
  x =  (float*) mxGetData(plhs[5]);
  cudaMemcpy(st, d_st, minSize * sizeof(int),   cudaMemcpyDeviceToHost);
  cudaMemcpy(id, d_id, minSize * sizeof(int),   cudaMemcpyDeviceToHost);
  cudaMemcpy(x,   d_x, minSize * sizeof(float), cudaMemcpyDeviceToHost);


  cudaFree(d_ftype);
  cudaFree(d_err);
  cudaFree(d_st);
  cudaFree(d_id);
  cudaFree(d_x);
  cudaFree(d_counter);
  cudaFree(d_Params);

  cudaFree(d_dout);

  mxGPUDestroyGPUArray(data);
  mxGPUDestroyGPUArray(dataraw);
  mxGPUDestroyGPUArray(WtW);
  mxGPUDestroyGPUArray(datarez);
  mxGPUDestroyGPUArray(W);
  mxGPUDestroyGPUArray(U);
  mxGPUDestroyGPUArray(dW);
  mxGPUDestroyGPUArray(dU);
  
}
