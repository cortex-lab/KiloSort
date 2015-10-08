%%
spk = cell(Nfilt,1);
for iNN = 1:Nfilt
   spk{iNN} = int32(rez.st3pos(rez.st3pos(:,2)==iNN, 1)); 
end

%% create templates
spktimes = [];
for iNN = 1:Nfilt
    eval(sprintf('spktimes.temp_%d = spk{iNN} + 21;', iNN-1));
end
%%
save('C:\\DATA\\Phy\\MariusTest\\spiketimes.mat','-struct','spktimes', '-v7.3');
%%
chanMap0 = chanMap(connected>1e-6);
for i = 1:length(chanMap0)
    chanMap0(i) = chanMap0(i) - sum(chanMap0(i) > chanMap(connected<1e-6));
end

[~, invchanMap0] = sort(chanMap0);

templates = zeros(Nchan, nt0, Nfilt, 'single');
for iNN = 1:Nfilt
   templates(chanMap0,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
end

%%
save('C:\\DATA\\Phy\\MariusTest\\templates.mat', 'templates', '-v7.3')