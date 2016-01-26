addpath('C:\CODE\GitHub\KiloSort')

ops.Nfilt               = 768 ; % 768 number of filters to use 
ops.Nfilt0              = 768 ; % 768 number of filters to use 
ops.nfullpasses         = 6; % 6, number of complete passes through data

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

ops.Th           = [4 12 10]; %[1 8 6]; % 4 12 10 6, 10, 8
ops.lam          = [10 30 30]; %[.1 30 30]; %10, 30, 30 %5, .... 0.1;
ops.nannealpasses = 4; %4
ops.momentum      = [20 400]; %[100 500];

ops.nNeighPC    = 12; % number of channnels to mask the PCs
ops.nNeig       = 16; % number of neighboring templates to retain projections of
%%
ops.ForceMaxRAMforDat   = Inf;

fidname{1}  = '20141202_all_es';
fidname{2}  = '20150924_1_e';
fidname{3}  = '20150601_all_s';
fidname{4}  = '20150924_1_GT';
fidname{5}  = '20150601_all_GT';
fidname{6}  = '20141202_all_GT';

for idset = 4
    ops.shuffle_clusters = 1;
    
    clearvars -except fidname ops idset  tClu tRes time_run
    
    root        = 'C:\DATA\Spikes';
    fname       = sprintf('set%d//%s.dat', idset, fidname{idset});
    fnameTW     = 'temp_wh.dat'; % (residual from RAM) of whitened data
    ops.chanMap = 'forPRBimecToWhisper.mat';
    
    clear loaded
    load_data_buff; % loads data into RAM + residual data on SSD
    %
    clear initialized
    run_reg_mu2; % iterate the template matching (non-overlapping extraction)
    %
    fullMPMU; % extracts final spike times (overlapping extraction)
    %
    if idset<=3
        testCode;
    end
end
% clear DATA
% plot_final_waveforms;
%%
