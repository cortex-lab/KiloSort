% mexGPUall. For these to complete succesfully, you need to configure the
% Matlab GPU library first (see README files for platform-specific
% information)
if ispc
    mexcuda -largeArrayDims mexMPmuFEAT.cu
    mexcuda -largeArrayDims mexMPregMU.cu
    mexcuda -largeArrayDims mexWtW2.cu
else
    mex -largeArrayDims mexMPmuFEAT.cu
    mex -largeArrayDims mexMPregMU.cu
    mex -largeArrayDims mexWtW2.cu
end


