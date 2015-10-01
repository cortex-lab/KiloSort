How to add CUDA paths to the systems path: append to /home/marius/.bashrc then logout/login:

export CUDA_HOME=/usr/local/cuda-7.5 
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64 
export PATH=${CUDA_HOME}/bin:${PATH} 

Then run mexGPUall.m