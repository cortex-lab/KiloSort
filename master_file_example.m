% default options are in parenthesis after the comment

addpath(genpath('C:\bin\GitHub\KiloSort')) % path to kilosort folder
addpath(genpath('C:\bin\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath('configFiles')
%% 
run('configFiles/StandardConfig.m')
tic; % start timer

if ops.GPU
     % initialize GPU (will erase any existing GPU arrays)
    gpuDevice(1);
end

if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
%
[rez, DATA, uproj] = preprocessData(ops);
rez = fitTemplates(ops, rez, DATA, uproj); 
rez = fullMPMU(ops, rez, DATA);% extracts final spike times (overlapping extraction)

% posthoc merge templates (under construction)
%     rez = merge_posthoc2(rez);

% save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);
%%
