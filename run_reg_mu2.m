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
%     Winit = Winit0;
    
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
    
    mu = 7 * ones(Nfilt, 1, 'single');
    
    nspikes = zeros(Nfilt, Nbatch);
    lam =  ones(Nfilt, 1, 'single');
 
    freqUpdate = 50;
    dWUtot= zeros(nt0, Nchan, Nfilt, 'single');
    iUpdate = 1:freqUpdate:Nbatch;
    
    i = 1;
    initialized = 1;

    npm = ones(Nfilt, 1);
    dbins = zeros(100, Nfilt);
    dsum = 0;
    miniorder = repmat(iperm, 1, ops.nfullpasses);
%     miniorder = repmat([1:Nbatch Nbatch:-1:1], 1, ops.nfullpasses/2);    
end


%%
% pmi = exp(-1./exp(linspace(log(ops.momentum(1)), log(ops.momentum(2)), Nbatch*ops.nannealpasses)));
pmi = exp(-1./linspace(ops.momentum(1), ops.momentum(2), Nbatch*ops.nannealpasses));

% pmi  = linspace(ops.momentum(1), ops.momentum(2), Nbatch*ops.nannealpasses);
Thi  = linspace(ops.Th(1),                 ops.Th(2), Nbatch*ops.nannealpasses);
if ops.lam(1)==0
    lami = linspace(ops.lam(1), ops.lam(2), Nbatch*ops.nannealpasses); 
else
    lami = exp(linspace(log(ops.lam(1)), log(ops.lam(2)), Nbatch*ops.nannealpasses));
end
 
gpuDevice(1);
if Nbatch_buff<Nbatch
    fid = fopen(fullfile(root, fnameTW), 'r');
end

st3 = [];

msg = [];
fprintf('Time %3.0fs. Optimizing templates ...\n', toc)
while (i<=Nbatch * ops.nfullpasses+1)
    if i<Nbatch*ops.nannealpasses
        Th      = Thi(i);
        lam(:)  = lami(i);
        pm      = pmi(i);
    end
    Params = double([NT Nfilt Th maxFR 10 Nchan Nrank]);    
    %
    if i>1 &&  ismember(rem(i,Nbatch), iUpdate) %&& i>Nbatch
        dWUtotCPU = gather(dWUtot);
        ntot = sum(nspikes,2);

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
                mu(k) = sum(sum(wu.*squeeze(dWUtotCPU(:,:,k))))/npm(k);
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
                
        rez.errall(ceil(i/freqUpdate))          = nanmean(delta);

        plot(sort(mu))
        axis tight
        drawnow
        
        
        %
        if ops.shuffle_clusters && i>Nbatch && rem(rem(i,Nbatch), 400)==1
            %
            uu = Nbatch * dbins/dsum;
            nhist = 1:1:100;
            
            [score, iY1, mu1, mu2, u1, u2]   = split_clust(uu, nhist);
            [d2d, iY, drez]                 = distance_betwxt(W, U, mu, sum(uu,1));
            
            [dsort, isort] = sort(drez, 'ascend');
            nswitch = find(score'<dsort(1:numel(score)), 1);
            if isempty(nswitch)
                nswitch = numel(score)+1;
            end
            nswitch = nswitch - 1;
%             disp(nswitch)
            
            for k = 1:nswitch
                
%                 mu0 = mu(iY1(k));
               mu(isort(k)) = mu1(k);
               mu(iY1(k))   = mu2(k);
               
               dbins(:, isort(k)) = u1(:, k) * dsum/Nbatch;
               dbins(:, iY1(k))   = u2(:, k) * dsum/Nbatch;

               W(:,isort(k),:) = W(:,iY1(k),:);
               U(:,isort(k),:) = U(:,iY1(k),:);

               ratio = sum(u1(:, k),1)/(sum(u1(:, k),1) + sum(u2(:, k),1));
               npm(isort(k)) = 1; %ratio * nsp(iY1(k));
               npm(iY1(k))   = 1; %(1-ratio) * nsp(iY1(k));
               
               dWUtot(:,:,isort(k)) = 0; %ratio * mu1(k)/mu0 * dWUtot(:,:,iY1(k));
               dWUtot(:,:,iY1(k))   = 0; %(1-ratio) * mu2(k)/mu0 * dWUtot(:,:,iY1(k));
            end
        end
    end
    
%     ibatch = iperm(rem(i-1,Nbatch)+1);
    ibatch = miniorder(i);
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
    
    % nonlinearity on the data
%     dataRAW = 8*(2./(1+exp(-dataRAW/4)) - 1);
    
    data 	= dataRAW * U(:,:); 
    
    U0 = gpuArray(U);
    utu = gpuArray.zeros(Nfilt, 'single');
    for irank = 1:Nrank
        utu = utu + (U0(:,:,irank)' * U0(:,:,irank));
    end
    UtU = logical(utu);
    
    [dWU, st, id, x,Cost] = mexMPregMU(Params,dataRAW,W(:,:),data,UtU,mu, lam * (20./mu).^2);
    nsp = histc(id, 0:1:size(W,2));
    nsp = nsp(1:Nfilt);
    nspikes(:, ibatch) = nsp;
    nsp = nsp(:);
%     dWU = dWU ./ permute(repmat(max(.1,  nsp(1:Nfilt)), 1, nt0, Nchan), [2 3 1]);
    
    x = min(max(1, round(x)), 100);
    dbins = .9975 * dbins;
    for j = 1:length(st)
        dbins(x(j), id(j)+1) = dbins(x(j), id(j)+1) + .0025;
    end
    dsum = .9975 * dsum +  .0025;

    fact = permute(repmat(pm.^nsp, 1, nt0, Nchan), [2 3 1]);
    
    dWUtot = fact .* dWUtot + (1-fact) .* gather(dWU);
    npm    = pm.^nsp .* npm    + (1-pm.^nsp) .* nsp;
    
    delta(ibatch) = sum(Cost)/1e6;
    
    if rem(i,100)==1
        nsort = sort(sum(nspikes,2), 'descend');
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Time %2.2f, batch %d/%d, mu %2.2f, err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d\n', ...
            toc, i,Nbatch* ops.nfullpasses,nanmedian(mu(:)), nanmean(delta), sum(nspikes(:)), ...
            nsort(min(size(W,2), 100)), nsort(min(size(W,2), 200)), ...
                nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
        fprintf(msg);        
    end
    i = i+1;
end
if Nbatch_buff<Nbatch
    fclose(fid);
end

% final templates computed here


