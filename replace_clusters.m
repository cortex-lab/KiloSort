function [W, U, mu, dWUtot, dbins, nswitch] = ...
    replace_clusters(dWUtot,W,U, mu, dbins, dsum, Nbatch, mergeT, splitT, Winit, Uinit, muinit)

uu = Nbatch * dbins/dsum;
nhist = 1:1:100;
nSpikes = sum(uu,1);

[score, iY1, mu1, mu2, u1, u2]   = split_clust(uu, nhist);
[~, iY, drez]                 = distance_betwxt(W, U, mu, nSpikes);

[dsort, isort] = sort(drez, 'ascend');

nmerged = sum(dsort<mergeT);
nsplit = sum(score>splitT);

% nswitch = find(score'<dsort(1:numel(score)), 1);
% if isempty(nswitch)
%     nswitch = numel(score)+1;
% end
% nswitch = nswitch - 1;
%             disp(nswitch)

freeInd = find(nSpikes<100 | mu'<10 | isnan(mu'));

for k = 1:nmerged
    % merge the two clusters
    iMerged = iY(isort(k));
    wt = [nSpikes(iMerged); nSpikes(isort(k))];
    wt = wt/sum(wt);
    mu(iMerged) = [mu(iMerged) mu(isort(k))] * wt;
    W(:,iMerged,:)  = W(:,iMerged,:) * wt(1) + W(:,isort(k),:) * wt(2);
    U(:,iMerged,:)  = U(:,iMerged,:) * wt(1) + U(:,isort(k),:) * wt(2);
    
    W(:,isort(k),:)  = 1e-10;
    U(:,isort(k),:)  = 1e-10;
end


for k = 1:min(nmerged+numel(freeInd), nsplit)
    if k<=numel(freeInd)
        inew= freeInd(k);
    else
%         isplit = k - numel(freeInd);
        inew = isort(k - numel(freeInd));
    end
    
    mu0 = mu(iY1(k));
    
    % split the bimodal cluster, overwrite merged cluster
    mu(inew)     = mu1(k);
    mu(iY1(k))   = mu2(k);
    
    dbins(:, inew)     = u1(:, k) * dsum/Nbatch;
    dbins(:, iY1(k))   = u2(:, k) * dsum/Nbatch;
    
    W(:,inew,:) = W(:,iY1(k),:);
    U(:,inew,:) = U(:,iY1(k),:);
    
    %     ratio = sum(u1(:, k),1)/(sum(u1(:, k),1) + sum(u2(:, k),1));
%     npm(inew)     = 20; %ratio * nsp(iY1(k));
%     npm(iY1(k))   = 20; %(1-ratio) * nsp(iY1(k));
    
    dWUtot(:,:,inew)     = 20 * mu1(k)/mu0 * dWUtot(:,:,iY1(k)); %/npm(iY1(k));
    dWUtot(:,:,iY1(k))   = 20 *mu2(k)/mu0 * dWUtot(:,:,iY1(k)); %/npm(iY1(k));
end

d2d                 = pairwise_dists(W, U, mu, Winit, Uinit, muinit);
dmatch              = min(d2d, [], 1);

inovel = find(dmatch(1:1000)>.4);
inovel = inovel(randperm(numel(inovel)));

i0 = Inf;

for k = 1+min(nmerged+numel(freeInd), nsplit):nmerged+numel(freeInd)
    % add new clusters
    i0 = i0 + 1;
    if i0>numel(inovel)
        break;
    end
    if k<=numel(freeInd)
        inew= freeInd(k);
    else
%         isplit = k - numel(freeInd);
        inew = isort(k - numel(freeInd));
    end
     
    mu(inew)           = muinit(inovel(i0));
    dbins(:, inew)     = 1;
    
    
    W(:,inew,:) = Winit(:,inovel(i0),:);
    U(:,inew,:) = Uinit(:,inovel(i0),:);
    
    %     ratio = sum(u1(:, k),1)/(sum(u1(:, k),1) + sum(u2(:, k),1));
    %npm(inew) = 20; %ratio * nsp(iY1(k));
    
    dWUtot(:,:,inew)     = 0; %ratio * mu1(k)/mu0 * dWUtot(:,:,iY1(k));
    for j = 1:size(W,3)
       dWUtot(:,:,inew) = dWUtot(:,:,inew) + W(:,inew,j) * U(:,inew,j)';
    end
    dWUtot(:,:,inew) = 20 * dWUtot(:,:,inew) * mu(inew);
end

nswitch = [min(nmerged, nsplit) i0]; %min(nmerged+numel(freeInd), nsplit);