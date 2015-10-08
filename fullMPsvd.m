% Params(3) = 6;
% Params(4) = 50000;
% Params(5) = 25; 

Params(2) = Nfilt;
Params(3) = 7;
Params(4) = 50000;
Params(5) = 10; 

U0 = gpuArray(U);
WtW  = gpuArray.zeros(Nfilt,Nfilt, 2*nt0-1, 'single');
for i = 1:Nrank
    for j = 1:Nrank
        utu0 = U0(:,:,i)' * U0(:,:,j);
        wtw0 = mexWtW2(Params, W(:,:,i), W(:,:,j), utu0);
%         wtw0 = squeeze(wtw(:,i,:,j,:));
        WtW = WtW + wtw0;
    end
end
mWtW = max(WtW, [], 3);
mWtW = gather(mWtW);

WtW = permute(WtW, [3 1 2]);

rez.WtW = gather(WtW);
clear wtw0 utu0 U0
%
clear nspikes2
st3 = [];

ipos = false(Nchan, Nfilt);
for i = 1:Nfilt
   ipos(:,i) = abs(U(:,i,1))>1e-3; 
   [~, sortU] = sort(abs(U(:,i,1)), 'descend');
   ipos(sortU(21:end),i) = false;
end

fid = fopen(fullfile(root, fnameTW), 'r');
msg = [];

Mask = gather(abs(WtW)>.01);
nt1 = nt0;
COV = cell(Nfilt,1);
for i = 1:Nfilt
   COV{i} = zeros(sum(ipos(:,i))*nt0, 'single'); 
end

if ops.verbose
   fprintf('Time %3.0fs. Running the final template matching pass...\n', toc) 
end
%%
for ibatch = 1:Nbatch
    if ibatch>Nbatch_buff
        offset = 2 * ops.Nchan*batchstart(ibatch-Nbatch_buff);
        fseek(fid, offset, 'bof');
        dat = fread(fid, [NT ops.Nchan], '*int16');
    else
       dat = DATA(:,:,ibatch); 
    end
    dataRAW = gpuArray(dat);
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    data 	= dataRAW * U(:,:); 
    
    [st, id, x] = mexMPmuLITE(Params,data,W(:,:),WtW, mu, lam * 20./mu);
    nspikes2(1:size(W,2)+1, ibatch) = histc(id, 0:1:size(W,2));

    inds = repmat(double(st'), nt0, 1) + repmat((1:nt0)', 1, numel(st));
    datS = reshape(dataRAW(inds, :), nt0, numel(st), Nchan);
    datS = permute(datS, [1 3 2]);
    
    for k = 1:Nfilt
        isk = (id==(k-1));
        if sum(isk)>0
            uu = datS(:, ipos(:,k), isk);
            uu = reshape(uu, [], size(uu,3))/50;
            COV{k} = COV{k} + gather(uu(:,:) * uu(:,:)');
        end
    end
    clear datS uu
    
    if rem(ibatch,100)==1
        nsort = sort(sum(nspikes2,2), 'descend');
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Time %2.2f, batch %d/%d, err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d\n', ...
            toc, ibatch,Nbatch, nanmean(delta), sum(nspikes2(:)), nsort(100), nsort(200), ...
            nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
        fprintf(msg);
    end
end

nsort = sort(sum(nspikes2,2), 'descend');
fprintf('Time %3.0fs. ExpVar %2.6f, n10 %d, n20 %d, n30 %d, n40 %d \n', toc, nanmean(delta), nsort(10), nsort(20), ...
    nsort(min(size(W,2), 30)), nsort(min(size(W,2), 40)));

fclose(fid);
%%
U3 = cell(Nfilt,1);
for i = 1:Nfilt
   [Usvd, Sv] = eigs(double(COV{i}), 10);
   U3{i} = Usvd(:,1:10);
end
%%

