ops.spkTh = -4;
ops.loc_range = [3  1];
ops.long_range = [30  6];
ops.maskMaxChannels = 5;

ops.crit = .1;
ops.nFiltMax = 10000;

dd = load('PCspikes.mat');
ops.wPCA = dd.Wi;

%%
tic

% find isolated spikes
[row, col, mu] = isolated_peaks(S1, loc_range, long_range, Th);

% find their PC projections
uS = get_PCproj(S1, row, col, wPCA, maskMaxChannels);

% merge in with existing templates
[nSnew, iNonMatch] = merge_spikes_in(uBase(:,1:ncurr,:), nS(1:ncurr), uS, crit);
nS(1:ncurr) = nS(1:ncurr) + nSnew;

% reduce non-matches
[uNew, nSadd] = reduce_clusters(uS(:,iNonMatch,:), crit);

% add new spikes to list
uBase(:, ncurr + [1:size(uNew,2)], :) = uNew;
nS(ncurr + [1:size(uNew,2)]) = nSadd;

ncurr = ncurr + size(uNew,2);

toc
%%

