function data = shift_data(data, dy, ycoords, xcoords, iCovChans, sigDrift, Wrot)

if nargin>6
    shiftM = shift_matrix(dy, ycoords, xcoords, iCovChans, sigDrift, Wrot);
else
    shiftM = shift_matrix(dy, ycoords, xcoords, iCovChans, sigDrift);
end

 data   = shiftM * data;



