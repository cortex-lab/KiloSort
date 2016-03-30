addpath('C:\CODE\GitHub\KiloSort')

ops.Nfilt               = 1024 ; %  number of filters to use (512)
ops.Nrank               = 3;    % matrix rank of spike template model
ops.nfullpasses         = 6;    % number of complete passes through data (6)

ops.whitening 	= 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
ops.maxFR       = 20000;  % maximum number of spikes to extract per batch (20000)

ops.fs          = 25000; % sampling rate


ops.fshigh      = 300;   % frequency for high pass filtering
ops.NchanTOT    = 129;   % total number of channels
ops.Nchan       = 120;   % number of active channels 
ops.ntbuff      = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.scaleproc   = 200;   % int16 scaling of whitened data
ops.verbose     = 1;     
ops.NT          = 32*1024+ ops.ntbuff;% this is the batch size, very important for memory reasons. 
% should be multiple of 32 (or higher power of 2) + ntbuff

% these options can improve/deteriorate results. when multiple values are 
% provided for an option, the first two are beginning and ending anneal values, 
% the third is the value used in the final pass. 
ops.Th               = [4 12 12];    % threshold for detecting spikes on template-filtered data ([6 12 12])
ops.lam              = [10 30 30];   % large means amplitudes are forced around the mean ([10 30 30])
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)
ops.momentum         = 1./[20 400]; % start with high momentum and anneal (1./[20 1000])
ops.shuffle_clusters = 1;            % allow merges and splits (1)
ops.mergeT           = .1;          % upper threshold for merging (.1)
ops.splitT           = .1;           % lower threshold for splitting (.1)

ops.nNeighPC    = 12; %12; % number of channnels to mask the PCs, leave empty to skip (12)
ops.nNeigh      = 16; % number of neighboring templates to retain projections of (16)

% new options
ops.initialize = 'fromData'; %'fromData'; %'fromData';

% options for initializing spikes from data
ops.spkTh           = -4;      % spike threshold in standard deviations (4)
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)
dd                  = load('PCspikes.mat'); % you might want to recompute this from your own data
ops.wPCA            = dd.Wi;   % PCs 

%%
addpath('D:\DATA\Spikes\EvaluationCode')

ops.ForceMaxRAMforDat   = Inf; %Inf;

fidname{1}  = '20141202_all_es';
fidname{2}  = '20150924_1_e';
fidname{3}  = '20150601_all_s';
fidname{4}  = '20150924_1_GT';
fidname{5}  = '20150601_all_GT';
fidname{6}  = '20141202_all_GT';

for idset = 6
    
    clearvars -except fidname ops idset  tClu tRes time_run
    
    root        = 'C:\DATA\Spikes';
    fname       = sprintf('set%d//%s.dat', idset, fidname{idset});
    fnameTW     = 'temp_wh.dat'; % (residual from RAM) of whitened data
    ops.chanMap = 'C:\DATA\Spikes\forPRBimecToWhisper.mat';
    
    
    clear loaded
    load_data_and_initialize; % loads data into RAM + residual data on SSD
    
   %%
    clear initialized
    
%     optimizePeaks;
    run_reg_mu2; % iterate the template matching (non-overlapping extraction)
    %
    fullMPMU; % extracts final spike times (overlapping extraction)
    %
    testCode;
end
% clear DATA
% plot_final_waveforms;
%%
