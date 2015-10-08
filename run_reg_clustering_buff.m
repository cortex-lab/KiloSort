if ~exist('initialized')
    addpath('C:\CODE\MariusBox\Primitives\')

    Nfilt 	= ops.Nfilt; %256+128;
    nt0 	= 61;
    ntbuff    = ops.ntbuff;
    NT  	= 128*1024+ ntbuff;

    Th 		= ops.Th;
    maxFR 	= ops.maxFR;
   
    Nchan 	= ops.Nchan;

    batchstart = 0:NT:NT*(Nbatch-Nbatch_buff);
    
    delta = NaN * ones(Nbatch, 1);
    iperm = randperm(Nbatch);

    Params = double([NT Nfilt 10 maxFR 10 Nchan]);
    
    initialize_waves0;
    ipck = randperm(size(Winit,2), Nfilt);
    W = Winit(:, ipck);
    U = Uinit(:, ipck);

    dW0 = zeros(size(W), 'single');
    dU0 = zeros(size(U), 'single');

   
    nspikes = zeros(Nfilt+1, Nbatch);
    lam = 0 * ops.lam * ones(Nfilt, 1, 'single');
 
    freqUpdate = 100;
    NbinsUpdate = ceil(Nbatch/freqUpdate);
    dWUtot= gpuArray.zeros(nt0, Nchan, Nfilt, NbinsUpdate, 'single');
    iUpdate = 1:freqUpdate:Nbatch;
    
    i = 0;
    initialized = 1;
end

%%
fid = fopen(fullfile(root, fnameTW), 'r');

msg = [];
fprintf('Time %3.0fs. Optimizing templates ...\n', toc)
while (i<Nbatch * ops.nfullpasses)
    i = i+1;
    if i>Nbatch && ismember(rem(i,Nbatch), iUpdate)
        dWUtotCPU = gather(sum(dWUtot, 4));
        for k = 1:Nfilt
            [Uall, Sv, Vall] = svd(gather(dWUtotCPU(:,:,k)), 0);
            
            [~, imax] = max(abs(Uall(:,1)), [], 1);
            
            W(:,k) = - Uall(:,1)* sign(Uall(imax,1));
            U(:,k) = - Vall(:,1) * sign(Uall(imax,1));
        end
        
        ibacurr = rem(i-1, Nbatch)+1;
        dWUtot(:,:,:,ceil(ibacurr/freqUpdate)) = 0;
        rez.errall(ceil(i/freqUpdate)) = nanmean(delta);
        %
        mmax = max(U,[],1);
        U(abs(U)<.1*repmat(mmax, Nchan,1)) = 0;
        U = normc(U);
    end
    
    ibatch = iperm(rem(i-1,Nbatch)+1);
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
    data 	= dataRAW * U; 
    
    U0 = gpuArray(U);
    utu = U0' * U0;
    WtW = mexWtW(Params, W, utu);
    WtW = permute(WtW, [3 1 2]);
  
    UtU = logical(utu);
    
    [dWU, st, id, x] = mexMPreg(Params,dataRAW,W,data, UtU);
    
    ibacurr = rem(i-1, Nbatch)+1;
    dWUtot(:,:,:,ceil(ibacurr/freqUpdate)) = ...
        dWUtot(:,:,:,ceil(ibacurr/freqUpdate)) + dWU/1e4;
    
    nspikes(1:size(W,2)+1, ibatch) = histc(id, 0:1:size(W,2));
    delta(ibatch) = sum(x.^2)/1e6;
    
    if rem(i,100)==1
        nsort = sort(sum(nspikes,2), 'descend');
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Time %2.2f, batch %d/%d, err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d', ...
            toc, i,Nbatch * ops.nfullpasses, nanmean(delta), sum(nspikes(:)), nsort(100), nsort(200), ...
                nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
        fprintf(msg);        
    end
    
end
fprintf(repmat('\b', 1, numel(msg)));
msg = sprintf('Time %2.2f, batch %d/%d, err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d', ...
    toc, i,Nbatch * ops.nfullpasses, nanmean(delta), sum(nspikes(:)), nsort(100), nsort(200), ...
    nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
fprintf(msg);        
fprintf('\n')

% final templates computed here
dWUtotCPU = gather(sum(dWUtot, 4));
for k = 1:Nfilt
    [Uall, Sv, Vall] = svd(gather(dWUtotCPU(:,:,k)), 0);
    
    [~, imax] = max(abs(Uall(:,1)), [], 1);
    
    W(:,k) = - Uall(:,1)* sign(Uall(imax,1));
    U(:,k) = - Vall(:,1) * sign(Uall(imax,1));
end


fclose(fid);

% max(x(1:end))
%numel(st(:))


%clear dout dout2

