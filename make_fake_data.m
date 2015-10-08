rezGT = rez;
rezGT.st3 = [];
rezGT.st3pos = [];
rezGT.WtW = [];
%% make_fake_traces
nGT = 5;

rezGT.Ugt = zeros(Nchan, nGT, Nrank);
rezGT.Wgt = [];
rezGT.Sgt = [];

rng(1);

iNN = find(rezGT.nbins>1000);
iperm = randperm(numel(iNN));


igather = 0;
for i = 1:numel(iNN)
    imin = find(abs(rezGT.U(:,iNN(i),1))>.025, 1,'first');
    imax = find(abs(rezGT.U(:,iNN(i),1))>.025, 1,'last');
    
    if imax-imin>20 || imax-imin<5
        continue;
    end
    
    igather = igather + 1;
    if (imax+imin)/2 < 60
        rezGT.Ugt((imin:imax) + 20 + ceil(rand*20),igather,:) = rezGT.U(imin:imax,iNN(i),:);
        rezGT.Wgt(:,igather,:) = rezGT.W(:,iNN(i),:);
    else
        rezGT.Ugt((imin:imax) - 20 - ceil(rand*20),igather,:) = rezGT.U(imin:imax,iNN(i),:);
        rezGT.Wgt(:,igather,:) = rezGT.W(:,iNN(i),:);
    end
    
    if igather==nGT
       break; 
    end
end

%
rng(1);
FR      = 10;
muisi   = ops.fs/FR;

for i = 1:nGT
    isi = 42 + exprnd(muisi - 42, 1, ceil(2*NT*Nbatch/muisi));
    spks = round(cumsum(isi));
    spks(spks>(NT-ops.ntbuff)*size(DATA,3) - 64) = [];
    rezGT.Sgt{i} = spks;    
end
%% add spikes
MUgt    = 20;

for i = 1:nGT    
    Ugt = squeeze(rezGT.Ugt(:, i, :));
    Wgt = squeeze(rezGT.Wgt(:, i, :));
    
    wu = int16(MUgt * Wgt * Ugt' * ops.scaleproc);
    rezGT.wu{i} = wu;
    
    ibatch  = ceil(rezGT.Sgt{i}/(NT-ops.ntbuff));
    spbatch = rem(rezGT.Sgt{i}-1, (NT-ops.ntbuff))+1;
    for j = 1:length(rezGT.Sgt{i})
        DATA(spbatch(j) + (1:1:nt0), :, ibatch(j)) = ...
            DATA(spbatch(j) + (1:1:nt0), :, ibatch(j)) + ...
            wu;
    end
end

flag_added  = 1;
%% undo spikes
for i = 1:nGT    
    Ugt = squeeze(rezGT.Ugt(:, i, :));
    Wgt = squeeze(rezGT.Wgt(:, i, :));
    
    wu = int16(MUgt * Wgt * Ugt' * ops.scaleproc);
    rezGT.wu{i} = wu;
    
    ibatch  = ceil(rezGT.Sgt{i}/(NT-ops.ntbuff));
    spbatch = rem(rezGT.Sgt{i}-1, (NT-ops.ntbuff))+1;
    for j = 1:length(rezGT.Sgt{i})
        DATA(spbatch(j) + (1:1:nt0), :, ibatch(j)) = ...
            DATA(spbatch(j) + (1:1:nt0), :, ibatch(j)) - ...
            wu;
    end
end

flag_added  = 0;

%%
