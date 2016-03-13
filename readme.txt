Assuming gpu functions have been correctly compiled (see below), a "master_file.m" is available that you should copy to a local path and change for each of your experiments. The logic is that the git folder might be updated, and when that happens all extraneous files in that folder will be deleted and any changes you made reverted. 

The following are instructions for setting up mex compilation of CUDA files with direct Matlab inputs. Note you need Matlab with the parallel computing toolbox. Instructions below are for Linux and  Windows. Mac instructions are in a separate file (I haven't tested them though). If successful, you should be able to run mexGPUall. 


Linux

Install CUDA (should ask for a compatible recent version of gcc, will install Nvidia drivers if necessary).

Append to /home/marius/.bashrc then logout/login:
export CUDA_HOME=/usr/local/cuda-7.5 
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64 
export PATH=${CUDA_HOME}/bin:${PATH} 

Make sure a compilation options file mex_CUDA_glnxa64.xml is available in the code folder. 


Windows

This is more complicated in Windows as a combination of several factors
CUDA requires visual studio express to compile
64bit tools are only available for free in SDK 7.1 as an addon for Visual Express 2010, not later editions
Microsoft no longer supports Visual Studio 2010 (SDK 7.1 will not work in Windows 10 at all). 
SDK 7.1 does not create the required vcvars64.bat
Several library files are simply not copied by Matlab from the SDK into its own private folder

Steps 1-3 (to be completed in this order).
Install Visual Studio Expres 2010 in its default install location from a third party location like http://microsoft-visual-cpp-express.soft32.com/. Note if you don’t install in default location, it puts an important part of the installation on C: anyway and you may have to change default paths for the mex options file and vcvars64.bat file. 
Install Windows SDK 7.1. Make sure it’s installed in the same folder as Visual Studio 2010 because they share an environment variable VS100COMNTOOLS that the mex file uses to locate them. 
Install CUDA (comes with compatible Nvidia drivers). If you get an error of not finding the GPU at the beginning of installation, you should try newer Nvidia drivers. Ignoring the warning appears to install CUDA without proper paths and then compilation does not work. 

Step4. SDK 7.1 does not create the required vcvars64.bat in \VC\bin\amd64. Create a text file into another location (such as on the Desktop), put this line only into it
CALL C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd /x64
and save it as vcvars64.bat. Now copy this file into \VC\bin\amd64. 

Step5. mex_CUDA_win64.xml is the options file for CUDA compilation and it should have been included with the code you are trying to run. Otherwise, you have to copy it from a Matlab location into your code folder:
matlabroot/toolbox/distcomp/gpu/extern/src/mex/win64
Step6. Ensure mex_CUDA_win64.xml passes an option to the CUDA compiler that it should use VC 2010. This is a line that starts COMPFLAGS="--cl-version 2010 …”.  This is NOT the default for the file provided at the Matlab installation location. Additionally to speed up compilation time, remove all compile options for achitectures other than the one you are targeting (look it up online or identify it in Matlab with the gpuDevice command).
Step7. Windows libraries that Matlab does not install are available in the SDK installation . A copy of these files should have been included with this code compressed in root\libs. These have been taken from a path like
C:\Program Files\Microsoft SDKs\Windows\v7.1\Lib
and should be put into a path like
C:\Program Files\MATLAB\R2014b\extern\lib\win64\microsoft
If the libs have not been included, copy them over from the SDK. The file names are capitalized in the SDK but this does not matter. The missing library files are on the line of the SDK that starts like LINKLIBS="/LIBPATH and they are called:  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib  ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib

Step8. If your video card is also driving your display, you need to disable the Windows watchdog that kills any process occupying the GPU for too long. 
start regedit,
navigate to HKEY_LOCAL_MACHINE\System\CurrentControlSet\Control\GraphicsDrivers
create new DWORD key called TdrLevel, set value to 0,
restart PC.



