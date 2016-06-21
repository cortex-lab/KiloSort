% mexGPUall. For these to complete succesfully, you need to configure the
% Matlab GPU library first (see README files for platform-specific
% information)

mex -largeArrayDims mexWtW.cu
mex -largeArrayDims mexMPreg.cu
mex -largeArrayDims mexMPsub.cu
mex -largeArrayDims mexMPmuFEAT.cu
mex -largeArrayDims mexMPregMU.cu
mex -largeArrayDims mexWtW2.cu



