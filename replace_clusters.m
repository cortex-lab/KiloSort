function [W, U, npm, dWUtot, dbins, nswitch] = ...
    replace_clusters(dWUtot,W,U, npm, mu, dbins, dsum, Nbatch, mergeT, splitT)

uu = Nbatch * dbins/dsum;
nhist = 1:1:100;
nSpikes = sum(uu,1);

[score, iY1, mu1, mu2, u1, u2]   = split_clust(uu, nhist);
[d2d, iY, drez]                 = distance_betwxt(W, U, mu, nSpikes);

[dsort, isort] = sort(drez, 'ascend');

nmerged = sum(dsort<mergeT);
nsplit = sum(score>splitT);

% nswitch = find(score'<dsort(1:numel(score)), 1);
% if isempty(nswitch)
%     nswitch = numel(score)+1;
% end
% nswitch = nswitch - 1;
%             disp(nswitch)

freeInd = find(nSpikes<100 | mu'<10);

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
    
    
    % split the bimodal cluster, overwrite merged cluster
    mu(inew) = mu1(k);
    mu(iY1(k))   = mu2(k);
    
    dbins(:, inew)     = u1(:, k) * dsum/Nbatch;
    dbins(:, iY1(k))   = u2(:, k) * dsum/Nbatch;
    
    W(:,inew,:) = W(:,iY1(k),:);
    U(:,inew,:) = U(:,iY1(k),:);
    
    %     ratio = sum(u1(:, k),1)/(sum(u1(:, k),1) + sum(u2(:, k),1));
    npm(inew) = 1; %ratio * nsp(iY1(k));
    npm(iY1(k))   = 1; %(1-ratio) * nsp(iY1(k));
    
    dWUtot(:,:,inew)     = 0; %ratio * mu1(k)/mu0 * dWUtot(:,:,iY1(k));
    dWUtot(:,:,iY1(k))   = 0; %(1-ratio) * mu2(k)/mu0 * dWUtot(:,:,iY1(k));
end

for k = 1+min(nmerged, nsplit):nmerged
    % add new clusters
end

nswitch = min(nmerged+numel(freeInd), nsplit);