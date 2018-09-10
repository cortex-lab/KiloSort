function WUinit = optimizePeaks(ops, uproj)

nt0             = ops.nt0;

nProj           = size(uproj,2);
nspike = size(uproj, 1);
nSpikesPerBatch = 4000;
nbatch = ceil(nspike / nSpikesPerBatch);


uBase = zeros(1e4, nProj);
nS = zeros(1e4, 1);
ncurr = 1;

for ibatch = 1:nbatch
    % merge in with existing templates
    [b, e] = range_idx(ibatch, nSpikesPerBatch, nspike);
    uS = uproj(b:e, :);
    [nSnew, iNonMatch] = merge_spikes0(uBase(1:ncurr,:), nS(1:ncurr), uS, ops.crit);
    nS(1:ncurr) = nSnew;

    % reduce non-matches
    [uNew, nSadd] = reduce_clusters0(uS(iNonMatch,:), ops.crit);

    % add new spikes to list
    uBase(ncurr + [1:size(uNew,1)], :) = uNew;
    nS(ncurr + [1:size(uNew,1)]) = nSadd;

    ncurr = ncurr + size(uNew,1);

    if ncurr>1e4
        break;
    end
end

nS = nS(1:ncurr);
uBase = uBase(1:ncurr, :);
[~, itsort] = sort(nS, 'descend');

%% initialize U
Nfilt = ops.Nfilt;
lam = ops.lam(1) * ones(Nfilt, 1, 'single');

ind_filt = itsort(rem([1:Nfilt]-1, numel(itsort)) + 1);
if ops.GPU
    U = gpuArray(uBase(ind_filt, :))';
else
    U = uBase(ind_filt, :)';
end
U = U + .001 * randn(size(U));
mu = sum(U.^2,1)'.^.5;
U = normc(U);

for i = 1:10

    dWU = zeros(Nfilt, nProj, 'single');
    if ops.GPU
        nToT = gpuArray.zeros(Nfilt, 1, 'single');
    else
        nToT = zeros(Nfilt, 1, 'single');
    end

    for ibatch = 1:nbatch
        % find clusters
        [b, e] = range_idx(ibatch, nSpikesPerBatch, nspike);
        arr_rows = e - b + 1;
        if ops.GPU
            clips = reshape(gpuArray(uproj(b:e, :)), arr_rows, nProj);
        else
            clips = reshape(uproj(b:e, :), arr_rows, nProj);
        end
        ci = clips * U;
        
        ci = bsxfun(@plus, ci, (mu .* lam)');
        cf = bsxfun(@rdivide, ci.^2, 1 + lam');
        cf = bsxfun(@minus, cf, (mu.^2.*lam)');

        [~, id] = max(cf, [], 2);

        id = gather_try(id);

        if ops.GPU
            L = gpuArray.zeros(Nfilt, arr_rows, 'single');
        else
            L = zeros(Nfilt, arr_rows, 'single');
        end
        L(id.' + (0:Nfilt:(Nfilt * arr_rows - 1))) = 1;
        dWU = dWU + L * clips;
        
        nToT = nToT + sum(L, 2);
    end
    dWU  = bsxfun(@rdivide, dWU, nToT);
    
    U = dWU';
    mu = sum(U.^2,1)'.^.5;
    U = normc(U);
end

Nchan = ops.Nchan;
Nfilt = ops.Nfilt;
Nrank = ops.Nrank;
wPCA = ops.wPCA(:,1:3);
Urec = reshape(U, Nchan, size(wPCA,2), Nfilt);

Urec= permute(Urec, [2 1 3]);
Wrec = reshape(wPCA * Urec(:,:), nt0, Nchan, Nfilt);

Wrec = gather_try(Wrec);

W = zeros(nt0, Nfilt, Nrank, 'single');
U = zeros(Nchan, Nfilt, Nrank, 'single');

Wrec(isnan(Wrec(:))) = 0;
for j = 1:Nfilt
    [w, sv, u] = svd(Wrec(:,:,j));
    w = w * sv;
    
    Sv = diag(sv);
    W(:,j,:) = w(:, 1:Nrank)/sum(Sv(1:ops.Nrank).^2).^.5;
    U(:,j,:) = u(:, 1:Nrank);
end

mu = gather_try(single(mu));
muinit = mu;

WUinit = zeros(nt0, Nchan, Nfilt);
for j = 1:Nfilt
    WUinit(:,:,j) = muinit(j)  * Wrec(:,:,j);
end
WUinit = single(WUinit);

end

function [b, e] = range_idx(idx, binsize, maxidx)
    b = min([(idx - 1) * binsize  + 1, maxidx]);
    e = min([idx * binsize, maxidx]);
end
