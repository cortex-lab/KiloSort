function [imax, amax] = max_interpolate3(amps, rez, sigDrift)

% sigDrift = 25;
% Wrot = rez.Wrot;
xcoords = rez.xc;
ycoords = rez.yc;

[amax, imax] = max(amps, [], 1);

[uniqy, ~, iy] = unique(ycoords);
[uniqx, ~, ix] = unique(xcoords);









