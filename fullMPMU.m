lam(:)    = ops.lam(3);
Params(3) = ops.Th(3);
Params(4) = 50000;
Params(5) = 50; 

U0 = gpuArray(U);
WtW  = zeros(Nfilt,Nfilt, 2*nt0-1, 'single');
for i = 1:Nrank
    for j = 1:Nrank
        utu0 = U0(:,:,i)' * U0(:,:,j);
        wtw0 =  gather(mexWtW2(Params, W(:,:,i), W(:,:,j), utu0));
        WtW = WtW + wtw0;
        clear wtw0 utu0
%         wtw0 = squeeze(wtw(:,i,:,j,:));
        
    end
end

WtW = gpuArray(WtW);
%%

mWtW = max(WtW, [], 3);
murep = repmat(mu, 1, Nfilt);
mWtW = mWtW .* (min(murep , murep')./max(murep , murep'));
mWtW = gather(mWtW);

mWtW = mWtW - diag(diag(mWtW));

WtW = permute(WtW, [3 1 2]);
% rez.WtW = gather(WtW);
%%
clear wtw0 utu0 U0
%
clear nspikes2
st3 = [];
rez.st3 = [];

if ops.verbose
   fprintf('Time %3.0fs. Running the final template matching pass...\n', toc) 
end

if Nbatch_buff<Nbatch
    fid = fopen(fullfile(root, fnameTW), 'r');
end
msg = [];

if ~isempty(ops.nNeighPC)
    nNeighPC  = ops.nNeighPC;
    nNeigh    = ops.nNeigh;
    load PCspikes
    
    rez.cProj = zeros(5e6, nNeigh, 'single');
    rez.cProjPC = zeros(5e6, 3*nNeighPC, 'single');

    % sort pairwise templates
    cr    = gather(squeeze( WtW(nt0, :,:)));
    cr(isnan(cr)) = 0;
    
    % save full similarity score
    rez.simScore = cr;
    
    [~, iNgsort] = sort(cr, 1, 'descend');
    rez.simScore = cr;
    maskTT = zeros(Nfilt, 'single');
    rez.iNeigh = iNgsort(1:nNeigh, :);
    for i = 1:Nfilt
        maskTT(rez.iNeigh(:,i),i) = 1;
    end
    
    % sort best channels
    [~, iNch] = sort(abs(U(:,:,1)), 1, 'descend');
    maskPC = zeros(Nchan, Nfilt, 'single');
    rez.iNeighPC = iNch(1:nNeighPC, :);
    for i = 1:Nfilt
        maskPC(rez.iNeighPC(:,i),i) = 1;
    end
    maskPC = repmat(maskPC, 3, 1);
end

irun = 0;
i1nt0 = int32([1:nt0])';
%
for ibatch = 1:Nbatch    
    if ibatch>Nbatch_buff
        offset = 2 * ops.Nchan*batchstart(ibatch-Nbatch_buff); % - ioffset;
        fseek(fid, offset, 'bof');
        dat = fread(fid, [NT ops.Nchan], '*int16');
    else
       dat = DATA(:,:,ibatch); 
    end
    dataRAW = gpuArray(dat);
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    
     % project data in low-dim space 
    data = gpuArray.zeros(NT, Nfilt, Nrank, 'single');
    for irank = 1:Nrank
        data(:,:,irank) 	= dataRAW * U(:,:,irank); 
    end
    data = reshape(data, NT, Nfilt*Nrank);
    
%     [st, id, x] = mexMPmuLITE(Params,data,W,WtW, mu, lam * 20./mu);
    [st, id, x, errC, proj] = mexMPmuFEAT(Params,data,W,WtW, mu, lam .* (20./mu).^2);    

    if ~isempty(ops.nNeighPC)
        % PCA coefficients
        inds = repmat(st', nt0, 1) + repmat(i1nt0, 1, numel(st));
        datSp = reshape(dataRAW(inds, :), [size(inds) Nchan]);
        coefs = reshape(Wi' * reshape(datSp, nt0, []), size(Wi,2), numel(st), Nchan);
        coefs = reshape(permute(coefs, [3 1 2]), [], numel(st));
        coefs = coefs .* maskPC(:, id+1);
        iCoefs = reshape(find(maskPC(:, id+1)>0), 3*nNeighPC, []);
        rez.cProjPC(irun + (1:numel(st)), :) = gather(coefs(iCoefs)');
        
        % template coefficients
        proj = maskTT(:, id+1) .* proj;
        iPP = reshape(find(maskTT(:, id+1)>0), nNeigh, []);
        rez.cProj(irun + (1:numel(st)), :) = proj(iPP)';
    
        % increment number of spikes
        irun = irun + numel(st);
    end
    
    if ibatch==1; ioffset = 0;
    else ioffset = ops.ntbuff;
    end
    st = st - ioffset;
    
%     nspikes2(1:size(W,2)+1, ibatch) = histc(id, 0:1:size(W,2));
    STT = cat(2, 20 + double(st) +(NT-ops.ntbuff)*(ibatch-1), ...
        double(id)+1, double(x), ibatch*ones(numel(x),1));
    st3 = cat(1, st3, STT);
    
    if rem(ibatch,100)==1
%         nsort = sort(sum(nspikes2,2), 'descend');
        fprintf(repmat('\b', 1, numel(msg)));
%         msg = sprintf('Time %2.2f, batch %d/%d, err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d\n', ...
%             toc, ibatch,Nbatch, nanmean(delta), sum(nspikes2(:)), nsort(min(size(W,2), 100)),nsort(min(size(W,2), 200)), ...
%             nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
        msg = sprintf('Time %2.2f, batch %d/%d, err %2.6f, NTOT %d\n', ...
            toc, ibatch,Nbatch, nanmean(delta), size(st3,1));        
        fprintf(msg);
        
    end
end

[~, isort]      = sort(st3(:,1), 'ascend');
st3             = st3(isort,:);

rez.st3         = st3;
if ~isempty(ops.nNeighPC)
    % re-sort coefficients for projections
    rez.cProj(irun+1:end, :)    = [];
    rez.cProjPC(irun+1:end, :)  = [];
    rez.cProjPC                 = reshape(rez.cProjPC, size(rez.cProjPC,1), [], 3);
    
    rez.cProj       = rez.cProj(isort, :);
    rez.cProjPC     = rez.cProjPC(isort, :,:);
    
    
    % re-index the template coefficients
    for ik = 1:Nfilt
        iSp = rez.st3(:,2)==ik;
        OneToN = 1:nNeigh;
        [~, isort] = sort(rez.iNeigh(:,ik), 'ascend');
        OneToN(isort) = OneToN;
        rez.cProj(iSp, :) = rez.cProj(iSp, OneToN);
        
        OneToN = 1:nNeighPC;
        [~, isort] = sort(rez.iNeighPC(:,ik), 'ascend');
        OneToN(isort) = OneToN;
        rez.cProjPC(iSp, :,:) = rez.cProjPC(iSp, OneToN, :);
    end
    
    rez.cProjPC = permute(rez.cProjPC, [1 3 2]);
end

%%
% nsort = sort(sum(nspikes2,2), 'descend');
% fprintf('Time %3.0fs. ExpVar %2.6f, n10 %d, n20 %d, n30 %d, n40 %d \n', toc, nanmean(delta), nsort(10), nsort(20), ...
%     nsort(min(size(W,2), 30)), nsort(min(size(W,2), 40)));

rez.ops      = ops;

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

rez.nbins = histc(rez.st3(:,2), .5:1:Nfilt+1);

[~, rez.ypos] = max(rez.U(:,:,1), [], 1);
if Nbatch_buff<Nbatch
    fclose(fid);
end

% center the templates
rez.W = cat(1, zeros(nt0 - 40, Nfilt, Nrank), rez.W);
rez.WrotInv = (rez.Wrot/200)^-1;

%
