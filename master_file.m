addpath('C:\CODE\GitHub\KiloSort')

ops.Nfilt               = 512 ; % 768 number of filters to use 
ops.Nfilt0              = 512 ; % 768 number of filters to use 
ops.nfullpasses         = 6; % 6, number of complete passes through data
ops.ForceMaxRAMforDat   = Inf; % if you want to force it to use a swap file on disk and less RAM, or no RAM at all (set to 0). 


% you shouldn't need to change these options
ops.whitening 	= 'full'; % type of whitening, only full for now
ops.learnWU     = 1;      % whether to learn templates
ops.maxFR       = 20000;  % maximum number of spikes to extract per batch
ops.lambda      = 1e3;    % not used
ops.Th0         = -7; % not currently used
ops.ntbuff      = 64;  % samples of symmetrical buffer for whitening and spike detection
ops.scaleproc   = 200; % int16 scaling of whitened data
ops.Nrank       = 3;

% these options you should know how to change
ops.fs          = 25000; % sampling rate
ops.fshigh      = 300; % high-pass filtering the data above this frequency
ops.NchanTOT    = 129; % total number of channels
ops.Nchan       = 120; % number of active channels 
ops.verbose     = 1;

% the options might need to be experimented with. 
% where there are more than one value for a parameter (Th, lam, momentum) these indicate start and end values during optimization + (possibly) the value during final spike extraction

ops.Th            = [4 12 10];  % where to set the thresold for spikes 
ops.lam           = [10 30 30]; % how much to force spikes in the same cluster to have the same amplitude (
ops.momentum      = [20 400];   % how many detected spikes to average over (start and end values)
ops.nannealpasses = 4;          % annealing stops after this many passes

% feature extraction options. These features are used in the Phy GUI. 
ops.nNeighPC     = 12; % number of channnels to mask the PCs
ops.nNeigh       = 16; % number of neighboring templates to retain projections of

% data paths
root        = 'C:\DATA\Spikes';
% fname   = 'set2\20150924_1_e.dat';
% fname   = '20141202_all_es.dat';
fname       = 'set3/20150601_all_s.dat';
fnameTW     = 'temp_wh.dat'; % will be created. residual of whitened data (not fitting in RAM).
% if chanMap is a string, it will load the channel configuration from a file 
% if chanMap is an array, that will indicate the channel order indexing into the raw data channel order
% This can also be used to drop inactive channels. 
ops.chanMap = 'forPRBimecToWhisper.mat'; % [1:32]

%% if you need to reload the data, clear variable 'loaded'
clear loaded
load_data_buff; % loads data into RAM + residual data on SSD

%% if you need to rerun the optimization, clear variable 'initialized'
clear initialized
run_reg_mu; % iterate the template matching (non-overlapping extraction)

%%
fullMPMU; % extracts final spike times (overlapping extraction)


