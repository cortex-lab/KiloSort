function WU = alignWU(WU)

[nt0 , Nchan, Nfilt] = size(WU);
[~, imin] = min(reshape(WU, nt0*Nchan, Nfilt), [], 1);

iMinChan = ceil(imin/nt0);


% imin = rem(imin-1, nt0) + 1;

% [~, imax] = min(W, [], 1);
% dmax = -(imin - 20);
% dmax = min(1, abs(dmax)) .* sign(dmax);
 
dmax = zeros(Nfilt, 1);
for i = 1:Nfilt
    wu = WU(:,iMinChan(i),i);
%     [~, imin] = min(diff(wu, 1));
    [~, imin] = min(wu);
    dmax(i) = - (imin- 20);
    
    wu = zeros(nt0, Nchan);
    if dmax(i)>0        
        wu((dmax(i) + 1):nt0, :) = WU(1:nt0-dmax(i),:, i);        
    else        
        wu(1:nt0+dmax(i), :) = WU((1-dmax(i)):nt0,:, i);        
    end
    WU(:, :,i) = wu;
end




