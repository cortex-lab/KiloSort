addpath('C:\CODE\GitHub\KiloSort')

ops.Nfilt               = 768 ; % 768 number of filters to use 
ops.Nfilt0              = 768 ; % 768 number of filters to use 
ops.nfullpasses         = 6; % 6, number of complete passes through data
ops.ForceMaxRAMforDat   = Inf;


ops.whitening 	= 'full'; % type of whitening, only full for now
ops.learnWU     = 1; % whether to learn templates
ops.maxFR       = 20000; % maximum number of spikes to extract per batch
ops.lambda      = 1e3; % not used
ops.Th0         = -7; % not currently used

ops.fs          = 25000; % sampling rate
ops.fshigh      = 300;
ops.NchanTOT    = 129; % total number of channels
ops.Nchan       = 120; % number of active channels 
ops.ntbuff      = 64;  % samples of symmetrical buffer for whitening and spike detection
ops.scaleproc   = 200; % int16 scaling of whitened data
ops.fprate      = 0.1; % estimated false positive rate (is it a spike or not?)
ops.verbose     = 1;
ops.Nrank       = 3;

ops.Th           = [6 10 8]; %6, 10, 8
ops.lam          = [10 20 20]; %5, .... 0.1;
ops.nannealpasses = 4; %4
ops.momentum      = [20 200];
% addpath('C:\CODE\MariusBox\FSkilosort')
%%
root        = 'C:\DATA\Spikes';
% fname   = 'set2\20150924_1_e.dat';
% fname   = '20141202_all_es.dat';
fname       = 'set3/20150601_all_s.dat';
fnameTW     = 'temp_wh.dat'; % (residual from RAM) of whitened data
ops.chanMap = 'forPRBimecToWhisper.mat';

clear loaded
load_data_buff; % loads data into RAM + residual data on SSD
%%
clear initialized
run_reg_mu; % iterate the template matching (non-overlapping extraction)
%%
fullMPMU; % extracts final spike times (overlapping extraction)
%
testCode;
% clear DATA
% plot_final_waveforms;
%%
