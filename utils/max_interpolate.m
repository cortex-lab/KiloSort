function [imax, amax] = max_interpolate(amps, rez, sigDrift)

% sigDrift = 25;
% Wrot = rez.Wrot;
xcoords = rez.xc;
ycoords = rez.yc;

[uniqy, ~, iy] = unique(rez.yc);

[~, indmax] = max(amps, [], 1);

indmax = iy(indmax);

chanDists= bsxfun(@minus, rez.yc, rez.yc').^2 + bsxfun(@minus, rez.xc, rez.xc').^2;
iCovChans = my_inv(exp(-chanDists/(2*sigDrift^2)), 1e-6);


dy = linspace(-20, 20, 200)';

imax = zeros(size(amps,2), 1);
amax = zeros(size(amps,2), 1);

for i = 1:numel(uniqy)
   ix = indmax==i; 
   
   yminusy = bsxfun(@minus, uniqy(i) + [dy], ycoords').^2 + ...
       bsxfun(@minus, repmat(mean(xcoords), numel(dy), 1), xcoords').^2;

   newSamp = exp(- yminusy/(2*sigDrift^2));

   shiftM = newSamp * iCovChans;
%    shiftM = (shiftM/Wrot);
   
   interp_amps = shiftM * amps(:, ix); 
   
   [km, ki] = max(interp_amps, [], 1);
   
   imax(ix) = dy(ki) + uniqy(i);
   amax(ix) = km;
end
%

