% Params(3) = 6;
% Params(4) = 50000;
% Params(5) = 25; 

Params(2) = Nfilt;
Params(3) = 5;
Params(4) = 50000;
Params(5) = 50; 

% ParamsW = Params;
% ParamsW(2)= Nrank*Nfilt;
% utu = gpuArray.ones(Nrank*Nfilt, 'single');
% wtw = mexWtW(ParamsW, W(:,:), utu);
% wtw = reshape(wtw, Nfilt, Nrank, Nfilt, Nrank, 2*nt0-1);

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
%
rez.WtW = gather(WtW);
clear wtw0 utu0 U0
%
clear nspikes2
st3 = [];

if ops.verbose
   fprintf('Time %3.0fs. Running the final template matching pass...\n', toc) 
end

fid = fopen(fullfile(root, fnameTW), 'r');
msg = [];

% Mask = abs(mWtW-diag(diag(mWtW)))>.05;
Mask = gather(abs(WtW)>.01);
nt1 = nt0;

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
    
    [~, isort] = sort(st, 'ascend');
    st0 = st(isort);
    x0  = x(isort);
    id0 = 1+ id(isort);
    isiso = get_isolated(st0, id0, Mask, nt1);
  
%     SORT THESE CORRECTLY BEFORE RUNNING AGAIN
    
    inds0 = repmat(double(st0' + (id0'-1)*NT), nt0, 1) + repmat((1:nt0)', 1, numel(st0));
    
    coefs = zeros(numel(st0), Nrank);
    for irank = 1:Nrank
        inds = inds0 + (irank-1)*Nfilt*NT;
        ww = reshape(data(inds), nt0, []);
        coefs(:,irank) = gather(sum(W(:,id0,irank) .* ww, 1));
    end
    
    %%
    STT = cat(2, double(st0) +(NT-ops.ntbuff)*(ibatch-1), double(id0)+1, double(x0), double(isiso), coefs);
    st3 = cat(1, st3, STT);
    
%     keyboard;
    
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
fprintf('Time %3.0fs. Thresholding spikes at false positive rate...\n', toc) 
st3pos = [];
fprate = ops.fprate;
Thx = zeros(Nfilt,1);
for idd = 1:1:Nfilt
    ix = find(st3(:,2)==idd);
    xs = st3(ix, 3);
    
    Mu = 10*ops.Th;
    Nbins = 1000;
    
    bbins = linspace(0, Mu, Nbins);
    hpos = cumsum(hist(Mu - xs(xs>0), bbins));
    hneg = cumsum(hist(Mu + xs(xs<0), bbins));
    
    ifirst = find(hneg./hpos > fprate, 1);
    if isempty(ifirst)
        ifirst = numel(bbins);
    end
    Thx(idd) = Mu - bbins(ifirst);
    
    st3pos = cat(1, st3pos, st3(ix(xs>Thx(idd)), :));
end

[~, isort] = sort(st3pos(:,1), 'ascend');
st3pos = st3pos(isort,:);

rez.st3      = st3; 
rez.st3pos   = st3pos; 
rez.ops      = ops;

% WUnorms = sum(sum(dWUtotCPU.^2, 2), 1).^.5;
% rez.template = gather(dWUtotCPU ./ repmat(WUnorms, nt0, Nchan, 1));

rez.W = W;
rez.U = U;
rez.t2p = [];
for i = 1:Nfilt
    wav0 = W(:,i,1);
    wav0 = my_conv(wav0', .5)';
   [~, itrough] = min(wav0);
    [~, t2p] = max(wav0(itrough:end));
    rez.t2p(i,1) = t2p;
    rez.t2p(i,2) = itrough;   
end

rez.nbins = histc(rez.st3pos(:,2), .5:1:Nfilt+1);

[~, rez.ypos] = max(rez.U(:,:,1), [], 1);

% estimateErrors;
%%
