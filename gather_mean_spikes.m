tic
NchanTOT = ops.NchanTOT;

root        = 'F:\DATA\Spikes';
fname       = fullfile(root, sprintf('MLong//%s.dat', fidname{idset}));

d = dir(fullfile(fname));
ops.sampsToRead = floor(d.bytes/NchanTOT/2);


NT          = 128*1024+ ops.ntbuff;
NTbuff      = NT + 4*ops.ntbuff;
Nbatch      = ceil(d.bytes/2/NchanTOT /(NT-ops.ntbuff));

% load data into patches, filter, compute covariance, write back to
% disk

fprintf('Time %3.0fs. Loading raw data... \n', toc);
fid = fopen(fname, 'r');
ibatch = 0;
Nchan = ops.Nchan;

Nchans = ops.Nchan;
ts = [1:1:141]';

clear stimes
for iNN = 1:size(rez.W,2)
    stimes{iNN} = rez.st3pos(rez.st3pos(:,2)==iNN,1);
end
% stimes = gtimes;

s = cell(1, numel(stimes));
Wraw = zeros(141, Nchans, numel(stimes));
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
        
    
    if ibatch==1; ioffset = 0;
    else ioffset = ops.ntbuff;
    end
    %
    for iNN = 1:numel(stimes)
        st = stimes{iNN} + ioffset - (NT-ops.ntbuff)*(ibatch-1) - 60;
        st(st<0) = [];
        st(st>NT-ops.ntbuff-141) = [];
        
        if ~isempty(st)
            inds = repmat(st', 141, 1) + repmat(ts, 1, numel(st));
            
            spks = reshape(dataRAW(inds, :), 141, numel(st), Nchans);
            
            spks = gather(spks);
            s{iNN} = cat(3, s{iNN}, permute(spks, [1 3 2]));
            Wraw(:,:,iNN) = Wraw(:,:,iNN) + ...
                gather(squeeze(sum(spks,2)));
        end
    end
    
end

for iNN = 1:numel(stimes)
     Wraw(:,:,iNN) = Wraw(:,:,iNN)/numel(stimes{iNN});
end
fprintf('Time %3.2f. Mean waveforms computed... \n', toc);






