if ~exist('loaded', 'var')
    tic    
    load('C:\DATA\Spikes\forPRBimecToWhisper.mat');
    chanMap = chanMap(connected>1e-6);
    
    batch_path = fullfile(root, 'batches');
    if ~exist(batch_path, 'dir')
        mkdir(batch_path);
    end
    NchanTOT = ops.NchanTOT;
    
    d = dir(fullfile(root, fname));
    ops.sampsToRead = floor(d.bytes/NchanTOT/2);
    
    dmem         = memory;
    memfree      = 4 * 2^30;
    memallocated = min(ops.ForceMaxRAMforDat, dmem.MemAvailableAllArrays) - memfree;
    memallocated = max(0, memallocated);
    nint16s      = memallocated/2;
    
    NT          = 128*1024+ ops.ntbuff;
    NTbuff      = NT + 2*ops.ntbuff;
    Nbatch      = ceil(d.bytes/2/NchanTOT /(NT-ops.ntbuff));
    Nbatch_buff = floor(nint16s/ops.Nchan /(NT-ops.ntbuff));
    Nbatch_buff = min(Nbatch_buff, Nbatch-1);
    
     % load data into patches, filter, compute covariance, write back to
    % disk
    [b1, a1] = butter(3, 300/ops.fs, 'high');
    
    fprintf('Time %3.0fs. Loading raw data... \n', toc);
    fid = fopen(fullfile(root, fname), 'r');
    ibatch = 0;
    Nchan = ops.Nchan;
    CC = gpuArray.zeros( Nchan,  Nchan, 'single');
    
    if ~exist('DATA', 'var')
        DATA = zeros(NT, ops.Nchan, Nbatch_buff, 'int16');
    end
    
    while 1
        ibatch = ibatch + 1;
            
        offset = max(0, 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - ops.ntbuff));
        if ibatch==1
            ioffset = 0;
        else
            ioffset = ops.ntbuff;
        end
        fseek(fid, offset, 'bof');
        buff = fread(fid, [NchanTOT NTbuff], '*int16');
        
        if isempty(buff)
            break;
        end
        nsampcurr = size(buff,2);
        if nsampcurr<NTbuff
            buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
        end
        dataRAW = gpuArray(buff);
        dataRAW = dataRAW';
        dataRAW = single(dataRAW);
        dataRAW = dataRAW(:, chanMap);
        
        datr = filter(b1, a1, dataRAW);
        datr = flipud(datr);
        datr = filter(b1, a1, datr);
        datr = flipud(datr);
        
        CC = CC + (datr' * datr)/NT;
        
        if ibatch<=Nbatch_buff
            DATA(:,:,ibatch) = gather(int16( datr(ioffset + (1:NT),:)));
        end        
    end
    CC = CC / ibatch;
    fclose(fid);
    fprintf('Time %3.0fs. Channel-whitening filters computed. \n', toc);

    fprintf('Time %3.0fs. Loading raw data and applying filters... \n', toc);
    
    [E, D] 	= svd(CC);
    eps 	= 1e-6;
    Wrot 	= E * diag(1./(diag(D) + eps).^.5) * E';
    Wrot = ops.scaleproc * Wrot;
    %%
    ibatch = 0;
    fid = fopen(fullfile(root, fname), 'r');
    fidW = fopen(fullfile(root, fnameTW), 'w');
    
    while 1
        ibatch = ibatch + 1;
        if ibatch<=Nbatch_buff
            datr = single(gpuArray(DATA(:,:,ibatch)));
        else
            offset = max(0, 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - ops.ntbuff));
            if ibatch==1
                ioffset = 0; 
            else
                ioffset = ops.ntbuff;
            end
            fseek(fid, offset, 'bof');
            
            buff = fread(fid, [NchanTOT NTbuff], '*int16');
            if isempty(buff)
                break;
            end
            nsampcurr = size(buff,2);
            if nsampcurr<NTbuff
                buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
            end
            
            dataRAW = gpuArray(buff);
            dataRAW = dataRAW';
            dataRAW = single(dataRAW);
            dataRAW = dataRAW(:, chanMap);
            
            datr = filter(b1, a1, dataRAW);
            datr = flipud(datr);
            datr = filter(b1, a1, datr);
            datr = flipud(datr);
            
            datr = datr(ioffset + (1:NT),:);
        end
        
        datr    = datr * Wrot;
        
         if ibatch<=Nbatch_buff
             DATA(:,:,ibatch) = gather(datr);
         else
            datcpu  = gather(int16(datr));
            fwrite(fidW, datcpu, 'int16');
         end
    end
    
    fclose(fidW);
    fclose(fid);
    fprintf('Time %3.2f. Whitened data written to disk... \n', toc);
    fprintf('Time %3.2f. Preprocessing complete!\n', toc);
    
    loaded = 1;
end





