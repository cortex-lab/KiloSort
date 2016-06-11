tic
% if ~isempty(ops.chanMap)
%     load(ops.chanMap);
%     chanMapConn = chanMap(connected>1e-6);
% else
%     chanMapConn = 1:ops.Nchan;
% % end
% NchanTOT = ops.NchanTOT;
% 
% d = dir(fullfile(root, fname));
% ops.sampsToRead = floor(d.bytes/NchanTOT/2);
% 
% 
% NT          = 128*1024+ ops.ntbuff;
% NTbuff      = NT + 4*ops.ntbuff;
% Nbatch      = ceil(d.bytes/2/NchanTOT /(NT-ops.ntbuff));

% load data into patches, filter, compute covariance, write back to
% disk

% root        = 'F:\DATA\Spikes\GT32';
% fname       = fullfile(root, sprintf('%s.dat', fidname));
% NchanTOT = 32;
% Nchan    = 32;

% fname = 'F:\DATA\Spikes\set11\20150601_all_GT91.dat';
fname = 'F:\DATA\Spikes\set13\20141202_all_GT245.dat';
% fname = 'F:\DATA\Spikes\set11\20150601_all_GT192.dat';
NchanTOT = 120;
Nchan = 120;

chanMapConn = 1:Nchan;
fprintf('Time %3.0fs. Loading raw data... \n', toc);
fid = fopen(fname, 'r');
ibatch = 0;

ts = [1:1:61]';

clear stimes
% gtClu = Clu;
% gtRes = Res ;
[iClu] = unique(gtClu);
for iNN = 1:numel(iClu)
    stimes{iNN} = double(gtRes(gtClu==iClu(iNN)));
end

% for iNN = 1:size(rez.W,2)
%     stimes{iNN} = rez.st3(rez.st3(:,2)==iNN,1);
% end

Wraw = zeros(61, Nchan, numel(stimes));
%
while 1
    ibatch = ibatch + 1;
    
    offset = max(0, 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
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
    
    offset = (NT-ops.ntbuff)*(ibatch-1)-64 -64+20;
    
    buff = buff(chanMapConn,:)';
    %
    for iNN = 1:numel(stimes)
        st = stimes{iNN} - offset;
        st(st<0) = [];
        st(st>NT-ops.ntbuff) = [];
        
        if ~isempty(st)
            inds = repmat(st', 61, 1) + repmat(ts, 1, numel(st));
            
            Wraw(:,:,iNN) = Wraw(:,:,iNN) + ...
                squeeze(sum(reshape(buff(inds, :), 61, numel(st), Nchan),2));
        end
    end
    
end
fclose(fid);

for iNN = 1:numel(stimes)
     Wraw(:,:,iNN) = Wraw(:,:,iNN)/numel(stimes{iNN});
end
fprintf('Time %3.2f. Mean waveforms computed... \n', toc);






