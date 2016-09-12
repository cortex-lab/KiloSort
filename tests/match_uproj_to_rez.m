function [iclust, pmin] = match_uproj_to_rez(rez, imax, tall)


[~, iTempMax] = max(sq(max(abs(rez.Wraw), [], 2)), [], 1);
% [~, iTempMax] = max(sq(max(abs(rez.U), [], 3)), [], 1);


% iTempMax(2:end) =iTempMax(1:end-1);
iTempMax = rez.yc(iTempMax);

tsort = tall;
imaxsort = imax;
% [tsort, isort] = sort(tall, 'ascend');
% imaxsort = imax(isort);

irez = 1;

nSpikes = size(rez.st3,1);
iclust = zeros(length(tsort), 1);
pmin   = zeros(length(tsort), 1);

for iproj = 1:length(tsort)
    while irez<=nSpikes && rez.st3(irez, 1)<tsort(iproj)
        irez = irez+1;
    end
    
    irange = irez + [-10:10];
    
    irange(irange<1 | irange>nSpikes) = [];
    
    pdist = 1 * abs(rez.st3(irange, 1) - tsort(iproj))/4 +  ...
        1 * abs(iTempMax(rez.st3(irange,2)) - imaxsort(iproj))/20;
    
    [pmin(iproj), imin] = min(pdist);
    
    iclust(iproj) = rez.st3(irange(imin), 2);
end



