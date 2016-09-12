function [imax, amax] = max_interpolate2(amps, rez, sigDrift)

% sigDrift = 25;
% Wrot = rez.Wrot;
xcoords = rez.xc;
ycoords = rez.yc;


[amax, imax] = max(amps, [], 1);

[uniqy, ~, iy] = unique(ycoords);
[uniqx, ~, ix] = unique(xcoords);

mask = exp( - bsxfun(@minus, ycoords, ycoords(imax)').^2/(2*sigDrift^2));

mask = mask .* amps; 
mask = bsxfun(@rdivide, mask, sum(mask,1));

imax = ycoords' * mask;







