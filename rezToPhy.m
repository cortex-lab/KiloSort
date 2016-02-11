
function [spikeTimes, clusterIDs, amplitudes, templates, templateFeatures, ...
    templateFeatureInds, pcFeatures, pcFeatureInds] = rezToPhy(rez, savePath)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)
%
% spikeTimes will be in samples, not seconds


spikeTimes = uint64(rez.st3(:,1));
% [spikeTimes, ii] = sort(spikeTimes);
clusterIDs = uint32(rez.st3(:,2));
amplitudes = rez.st3(:,3);


load(rez.ops.chanMap);
% chanMap0 = chanMap(connected>1e-6);
Nchan = rez.ops.Nchan;
nt0 = size(rez.W,1);
U = rez.U;
W = rez.W;

% for i = 1:length(chanMap0)
%     chanMap0(i) = chanMap0(i) - sum(chanMap0(i) > chanMap(connected<1e-6));
% end
% [~, invchanMap0] = sort(chanMap0);

templates = zeros(Nchan, nt0, rez.ops.Nfilt, 'single');
for iNN = 1:rez.ops.Nfilt
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatures = rez.cProj;
templateFeatureInds = uint32(rez.iNeigh);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

if ~isempty(savePath)
    
    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
    writeNPY(clusterIDs-1, fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
    writeNPY(templates, fullfile(savePath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
    Fs = rez.ops.fs;
    conn = logical(connected);
    chanMap0ind = int32(chanMap0ind);
    
    writeNPY(chanMap0ind(conn), fullfile(savePath, 'channel_map.npy'));
    %writeNPY(connected, fullfile(savePath, 'connected.npy'));
%     writeNPY(Fs, fullfile(savePath, 'Fs.npy'));
    writeNPY([xcoords(conn) ycoords(conn)], fullfile(savePath, 'channel_positions.npy'));
    
    writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
    writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
    writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
    writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    whiteningMatrix = rez.Wrot;
    writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
    
end
