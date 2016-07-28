% mexGPUall. For these to complete succesfully, you need to configure the
% Matlab GPU library first (see README files for platform-specific
% information)
if ispc
    mexcuda -largeArrayDims mexWtW.cu
    mexcuda -largeArrayDims mexMPreg.cu
    mexcuda -largeArrayDims mexMPsub.cu
    mexcuda -largeArrayDims mexMPmuFEAT.cu
    mexcuda -largeArrayDims mexMPregMU.cu
    mexcuda -largeArrayDims mexWtW2.cu
else
    mex -largeArrayDims mexWtW.cu
    mex -largeArrayDims mexMPreg.cu
    mex -largeArrayDims mexMPsub.cu
    mex -largeArrayDims mexMPmuFEAT.cu
    mex -largeArrayDims mexMPregMU.cu
    mex -largeArrayDims mexWtW2.cu
end


