Assuming gpu functions have been correctly compiled (see below), a "master_file.m" is available that you should copy to a local path and change for each of your experiments. The logic is that the git folder might be updated, and when that happens all extraneous files in that folder will be deleted and any changes you made reverted. 

The following are instructions for setting up mex compilation of CUDA files with direct Matlab inputs. Note you need Matlab with the parallel computing toolbox. Instructions below are for Linux and  Windows. Mac instructions are in a separate file (I haven't tested them though). If successful, you should be able to run mexGPUall. 

You should be able to run the code on a CPU without compiling any files, but it will be much, much slower than even on 250$ GPUs. For 32 channel data though, it might be fast enough. 

Windows

Install Visual Studio Community (2012 or 2013)
Install CUDA (comes with compatible Nvidia drivers, though it must be installed separately. If you need to download it, look here: https://developer.nvidia.com/cuda-downloads). If you get an error of not finding the GPU at the beginning of installation, you should try newer Nvidia drivers. Ignoring the warning appears to install CUDA without proper paths and then compilation does not work. 

Copy mex_CUDA_win64.xml (or nvcc_msvc120.xml, or a similarly named file, compatible with your Visual Studio installation version; 11 for 2012 and 12 for 2013) from here
matlabroot/toolbox/distcomp/gpu/extern/src/mex/win64
and into the KiloSort folder (or somewhere in your path). The included file with KiloSort will NOT be compatible with your local environment (unless you start changing paths inside it). 

If your video card is also driving your display, you need to disable the Windows watchdog that kills any process occupying the GPU for too long. 
start regedit,
navigate to HKEY_LOCAL_MACHINE\System\CurrentControlSet\Control\GraphicsDrivers
create new DWORD key called TdrLevel, set value to 0,
restart PC.

Linux

Install CUDA (should ask for a compatible recent version of gcc, will install Nvidia drivers if necessary).

Append to /home/marius/.bashrc then logout/login:
export CUDA_HOME=/usr/local/cuda-7.5 
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64 
export PATH=${CUDA_HOME}/bin:${PATH} 

Copy mex_CUDA_glnxa64.xml (or a similarly named file, compatible with your Visual Studio installation) from here
matlabroot/toolbox/distcomp/gpu/extern/src/mex/
and into the KiloSort folder (or somewhere in your path). The included file with KiloSort will NOT be compatible with your local environment (unless you start changing paths inside it). 


