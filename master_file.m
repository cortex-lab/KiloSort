
ops.Nfilt               = 512; % number of filters to use 
ops.nfullpasses         = 6; % number of complete passes through data
ops.ForceMaxRAMforDat   = Inf;


ops.whitening 	= 'full'; % type of whitening, only full for now
ops.learnWU     = 1; % whether to learn templates
ops.maxFR       = 20000; % maximum number of spikes to extract per batch
ops.Th          = 10; % spike detection threshold
ops.lambda      = 1e3; % not used
ops.Th0         = -7; % not currently used

ops.fs          = 25000; % sampling rate
ops.NchanTOT    = 129; % total number of channels
ops.Nchan       = 120; % number of active channels 
ops.ntbuff      = 64;  % samples of symmetrical buffer for whitening and spike detection
ops.lam         = 1;   % not used yet
ops.scaleproc   = 200; % int16 scaling of whitened data
ops.fprate      = 0.1; % estimated false positive rate (is it a spike or not?)
ops.verbose     = 1;

% addpath('C:\CODE\MariusBox\FSkilosort')
%%
root    = 'C:\DATA\Spikes';
fname   = '20150924_1.dat';
fnameTW = '20150924_1_tw.dat'; % (residual from RAM) of whitened data

load_data_buff; % loads data into RAM + residual data on SSD

run_reg_clustering_buff; % iterate the template matching (non-overlapping extraction)
 
fullMPbuff; % extracts final spike times (overlapping extraction)
%%
