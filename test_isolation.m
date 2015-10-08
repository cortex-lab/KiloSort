[~, isort] = sort(st, 'ascend');
st0 = st(isort);
x0= x(isort);
id0 = 1+ id(isort);
%%
mWtW0 = gather(mWtW);
mWtW0 = mWtW0 - diag(diag(mWtW0));
UtU = gather(UtU);
UtU = UtU - diag(diag(UtU));

nt1 = 41;

imin = 1;
icurrent = 1;
imax = 1;
nspbatch = numel(st0);

isol = true(nspbatch, 1);

tic
while(icurrent<=nspbatch)
    while (st0(imax) - st0(icurrent)) < nt1 && imax<nspbatch
        imax = imax+1;
    end
    
    while (st0(icurrent) - st0(imin)) > nt1
        imin = imin+1;
    end
    if sum(abs(mWtW0(id0(icurrent), id0(imin:imax-1)))>.01)
        isol(icurrent) = false;
    end
    
    icurrent = icurrent + 1;
end
toc
%%
iNN = iNN+1;

itrueISO = find(isol);
tsiso = double(st0(itrueISO(iNN)));
subplot(1,2,1)
imagesc(dataRAW(tsiso + (1:nt0), :)')
subplot(1,2,2)
k = id0(itrueISO(iNN));
wu = squeeze(W(:,k,:)) * squeeze(U(:,k,:))';
imagesc(wu')