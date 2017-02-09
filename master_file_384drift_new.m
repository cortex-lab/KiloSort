clear DATA

addpath(genpath('D:\CODE\GitHub\KiloSort')) % path to kilosort folder
addpath(genpath('C:\CODE\GitHub\npy-matlab')) % path to npy-matlab scripts

pathToYourConfigFile = 'D:\CODE\MariusBox\KiloSortLocal'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'configFileBench384.m'))

ops.fbinary  = 'F:\DATA\Spikes\Eijkman\2016-05-21\Eijkman_20160521_M2_g0_t0.imec_AP_CAR.bin';
ops.root     = fileparts(ops.fbinary);
ops.doDriftCorrection = 1;

tic; % start timer
if ops.GPU
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

if strcmp(ops.datatype , 'openEphys')
    ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
rez.ops = ops;
rez.ops.ForceMaxRAMforDat = 0;
rez.ops.Drift.tSmooth     = .25;
rez.ops.Drift.chSmooth    = 10;
rez.ops.spkTh             = -6;
rez.ops.initialize        = 'fromDriftCorrection';

if rez.ops.doDriftCorrection
    [rez, uprojDrift, indBatch] = collectRawClips(rez);
    rez = clusterAndDriftCorrection(rez, uprojDrift, indBatch);
end
clear uprojDrift

%
[rez, DATA, uproj] = preprocessData(rez); % preprocess data and extract spikes for initialization

%
% rez.ops.Nfilt = 512;
% rez.ops.initialize = 'no';

rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
%
rez.ops.nNeighPC   = 12;
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


rezToPhy(rez, ops.root);


    rez = merge_posthoc2(rez);

% save matlab results file
save(fullfile(ops.root,  sprintf('rez_new%d.mat', rez.ops.Nfilt)), 'rez', '-v7.3');
% remove temporary file
%     delete(ops.fproc);
