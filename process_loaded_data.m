if ~exist('loaded', 'var')
    tic 

    NchanTOT = ops.NchanTOT;
    
    NT          = size(DATA,1);
    NTbuff      = NT + 4*ops.ntbuff;
    Nbatch      = size(DATA,3);
    Nbatch_buff = Nbatch;
     %% load data into patches, filter, compute covariance, write back to
    % disk
    [b1, a1] = butter(3, ops.fshigh/ops.fs, 'high');
    
    Nchan = ops.Nchan;
    CC = gpuArray.zeros( Nchan,  Nchan, 'single');
    %%
    for ibatch = 1:Nbatch
        dataRAW = gpuArray(DATA(:,:,ibatch));
        dataRAW = single(dataRAW);

        datr = filter(b1, a1, dataRAW);
        datr = flipud(datr);
        datr = filter(b1, a1, datr);
        datr = flipud(datr);
                
        CC = CC + (datr' * datr)/NT;
        
       
        DATA(:,:,ibatch) = gather(int16( datr));
       
    end
    CC = CC / Nbatch;
    
    fprintf('Time %3.0fs. Channel-whitening filters computed. \n', toc);

    fprintf('Time %3.0fs. Loading raw data and applying filters... \n', toc);
    
    %%
    switch ops.whitening
        case 'diag'
            CC = diag(diag(CC));
    end
    
    [E, D] 	= svd(CC);
    eps 	= 1e-6;
    Wrot 	= E * diag(1./(diag(D) + eps).^.5) * E';
    Wrot    = ops.scaleproc * Wrot;
    %
    for ibatch=1:Nbatch
        datr = single(gpuArray(DATA(:,:,ibatch)));
        datr    = datr * Wrot;
        DATA(:,:,ibatch) = gather(int16(datr));
    end
    
    Wrot        = gather(Wrot);
    rez.Wrot    = Wrot;

    fprintf('Time %3.2f. Whitened data written to disk... \n', toc);
    fprintf('Time %3.2f. Preprocessing complete!\n', toc);
    
    loaded = 1;
end





