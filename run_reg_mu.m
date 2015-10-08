if ~exist('initialized')
    addpath('C:\CODE\MariusBox\Primitives\')
    rng(1);
    
    Nfilt 	= ops.Nfilt; %256+128;
    nt0 	= 61;
    ntbuff  = ops.ntbuff;
    NT  	= 128*1024+ ntbuff;

    Nrank   = ops.Nrank;
    Th 		= ops.Th;
    maxFR 	= ops.maxFR;
   
    Nchan 	= ops.Nchan;

    batchstart = 0:NT:NT*(Nbatch-Nbatch_buff);
    
    delta = NaN * ones(Nbatch, 1);
    iperm = randperm(Nbatch);
    
    initialize_waves0;
    ipck = randperm(size(Winit,2), Nfilt);
    W = [];
    U = [];
    for i = 1:Nrank
        W = cat(3, W, Winit(:, ipck)/Nrank);
        U = cat(3, U, Uinit(:, ipck));
    end
    W = alignW(W);
    for k = 1:Nfilt
            wu = squeeze(W(:,k,:)) * squeeze(U(:,k,:))';
            newnorm = sum(wu(:).^2).^.5;
            W(:,k,:) = W(:,k,:)/newnorm;
    end
    
    mu = 20 * ones(Nfilt, 1, 'single');
    
    nspikes = zeros(Nfilt+1, Nbatch);
    lam = 1 * ops.lam * ones(Nfilt, 1, 'single');
 
    freqUpdate = 100;
    NbinsUpdate = ceil(Nbatch/freqUpdate);

    dWUtot= zeros(nt0, Nchan, Nfilt, NbinsUpdate, 'single');
    mutot = zeros(Nfilt, NbinsUpdate, 'single');
    iUpdate = 1:freqUpdate:Nbatch;
    
    i = 1;
    initialized = 1;
end

%%
tic
Thi  = linspace(ops.Th0, ops.Th, Nbatch*ops.nannealpasses);
lami = exp(linspace(log(ops.lam0), log(ops.lam), Nbatch*ops.nannealpasses));

gpuDevice(1);
fid = fopen(fullfile(root, fnameTW), 'r');

msg = [];
fprintf('Time %3.0fs. Optimizing templates ...\n', toc)
while (i<=Nbatch * ops.nfullpasses+1)
    if i<Nbatch*ops.nannealpasses
        Th = Thi(i);
        lam(:) = lami(i);
    end
    Params = double([NT Nfilt Th maxFR 10 Nchan Nrank]);    
    
    if i>1 && ismember(rem(i,Nbatch), iUpdate) %&& i>Nbatch
        dWUtotCPU = gather(sum(dWUtot, 4));
        ntot = sum(nspikes,2);
        %         mu = gather(sum(mutot,2) ./ sum(nspikes(1:Nfilt,:),2));
        for k = 1:Nfilt
            if ntot(k)>5
                [Uall, Sv, Vall] = svd(gather(dWUtotCPU(:,:,k)), 0);
                
                Sv = diag(Sv);
                sumSv2 = sum(Sv(1:Nrank).^2).^.5;
                for irank = 1:Nrank
                    [~, imax] = max(abs(Uall(:,irank)), [], 1);
                    W(:,k,irank) = - Uall(:,irank) * sign(Uall(imax,irank)) * Sv(irank)/sumSv2;
                    U(:,k,irank) = - Vall(:,irank) * sign(Uall(imax,irank));
                end
                mmax = max(abs(U(:,k,1)));
                Usize = squeeze(abs(U(:,k,:)));
                Usize = Usize .* repmat(Sv(1:Nrank)'/Sv(1), Nchan, 1);
                ibad = max(Usize, [], 2) < .1 * mmax;
                
                U(ibad,k,:) = 0;
            end
        end
        
        for k = 1:Nfilt
            if ntot(k)>5                
                wu = squeeze(W(:,k,:)) * squeeze(U(:,k,:))';
                mu(k) = sum(sum(wu.*squeeze(dWUtotCPU(:,:,k))));
                mu(k) = mu(k)./(ntot(k)/1e4);
            end
        end
        if i<Nbatch * ops.nfullpasses
            W = alignW(W);
        end
        for k = 1:Nfilt
            if ntot(k)>5
                wu = squeeze(W(:,k,:)) * squeeze(U(:,k,:))';
                newnorm = sum(wu(:).^2).^.5;
                W(:,k,:) = W(:,k,:)/newnorm;
            end
        end
        
        if i>Nbatch * ops.nfullpasses
            break;
        end
        
        ibacurr = rem(i-1, Nbatch)+1;
        dWUtot(:,:,:,ceil(ibacurr/freqUpdate))  = 0;
%         mutot(:, ceil(ibacurr/freqUpdate))      = 0;
        rez.errall(ceil(i/freqUpdate))          = nanmean(delta);

        plot(sort(mu))
        axis tight
        drawnow
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
    data 	= dataRAW * U(:,:); 
    
    U0 = gpuArray(U);
    utu = gpuArray.zeros(Nfilt, 'single');
    for irank = 1:Nrank
        utu = utu + (U0(:,:,irank)' * U0(:,:,irank));
    end
    UtU = logical(utu);
    
    [dWU, st, id, x,Cost] = mexMPregMU(Params,dataRAW,W(:,:),data,UtU,mu, lam * 20./mu);
    
    ibacurr = rem(i-1, Nbatch)+1;
    dWUtot(:,:,:,ceil(ibacurr/freqUpdate)) = ...
        dWUtot(:,:,:,ceil(ibacurr/freqUpdate)) + gather(dWU/1e4);
    
%     mutot(:,ceil(ibacurr/freqUpdate)) =...
%         mutot(:,ceil(ibacurr/freqUpdate)) + gather(xhat);
    
    nspikes(1:size(W,2)+1, ibatch) = histc(id, 0:1:size(W,2));
    delta(ibatch) = sum(Cost)/1e6;
    
    if rem(i,100)==1
        nsort = sort(sum(nspikes,2), 'descend');
%         fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Time %2.2f, batch %d/%d, mu %2.2f, err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d\n', ...
            toc, i,Nbatch* ops.nfullpasses,nanmedian(mu(:)), nanmean(delta), sum(nspikes(:)), nsort(100), nsort(200), ...
                nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
        fprintf(msg);        
    end
    i = i+1;
end
fclose(fid);

% final templates computed here


fprintf(repmat('\b', 1, numel(msg)));
msg = sprintf('Time %2.2f, batch %d/%d, mu %2.2f, err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d\n', ...
    toc, i,Nbatch* ops.nfullpasses,nanmedian(mu(:)), nanmean(delta), sum(nspikes(:)), nsort(100), nsort(200), ...
    nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
fprintf(msg);
fprintf('\n')

