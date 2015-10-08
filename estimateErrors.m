%%
nfb = 32;
ParamsW     = Params;
ParamsW(2)  = Nfilt + nfb;

Wall        = cat(2, W, rezGT.Wgt);
Wall(:,Nfilt+32,:) = 0;
U0 = gpuArray(cat(2, U, rezGT.Ugt));
U0(:,Nfilt+32,:) = 0;

WtW  = gpuArray.zeros(Nfilt+nfb,Nfilt+nfb, 2*nt0-1, 'single');
for i = 1:Nrank
    for j = 1:Nrank
        utu0 = U0(:,:,i)' * U0(:,:,j);
        wtw0 = mexWtW2(ParamsW, Wall(:,:,i), Wall(:,:,j), utu0);
%         wtw0 = squeeze(wtw(:,i,:,j,:));
        WtW = WtW + wtw0;
    end
end

mWtW = squeeze(max(WtW, [], 3));
muall = cat(1, mu, MUgt * ones(5,1));
muall(Nfilt+nfb) = 0;

murep = repmat(muall, 1, Nfilt+nfb);
mWtW = mWtW .* (min(murep , murep')./max(murep , murep'));
mWtW = gather(mWtW);
mWtW = mWtW - diag(diag(mWtW));
clear wtw0 utu0 U0 WtW

itarget = Nfilt + (1:1:nGT);
%%
for iGT = 1:nGT
    mww             = mWtW(1:Nfilt, itarget(iGT));
    mww(isnan(mww)) = 0;
    [CCsort, isort] = sort(mww, 'descend');
    dmin            = Inf * ones(length(rezGT.Sgt{iGT}), 1);
    
    ikl = 1;
    tots = 0;
    
    if CCsort(1)>.85
       imatch =  isort(CCsort>.85);
    else
       imatch =  isort(1); 
    end
    
    ispk = ismember(rez.st3pos(:,2),imatch);
    
%     ispk = ispk & abs(rez.st3pos(:,3)-20)<2;
    
    spktimes    = rez.st3pos(ispk,1);
    fprintf('Neu %d Tmatch %d CC %2.2f Nmatches %d Nspikes %d Mu %2.2f \n', ...
        iGT, sum(CCsort>.9), CCsort(ikl), sum(CCsort>.95), numel(spktimes), mu(isort(ikl)));
    
    for i = 1:length(rezGT.Sgt{iGT})
        dmin(i) = min(dmin(i), min(abs(rezGT.Sgt{iGT}(i) - spktimes)));
    end
    ikl         = ikl + 1;
    tots = tots + numel(spktimes);
    
    
    det0(1,iGT) = mean(dmin<20);
    det0(2,iGT) = (tots - sum(dmin<20))/numel(rezGT.Sgt{iGT});
end
det0