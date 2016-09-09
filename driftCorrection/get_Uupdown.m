function [Uup, Udown] = get_Uupdown(U, sigShift, ops, rez, iCovChans, sigDrift)
% resample U with a shift up and a shift down (sigShift in um)
 
Uup   = shift_data(reshape(U, ops.Nchan, []),  sigShift, rez.yc, rez.xc, iCovChans, sigDrift, rez.Wrot);
Uup   = reshape(Uup, size(U));

Udown = shift_data(reshape(U, ops.Nchan, []), -sigShift, rez.yc, rez.xc, iCovChans, sigDrift, rez.Wrot);
Udown = reshape(Udown, size(U));