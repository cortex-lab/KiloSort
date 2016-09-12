% uprojDrift = reshape(uprojDrift, [], 120, 3);
% amps = log(1+sum(uprojDrift.^2,3))';
% tall = tsDrift;


uproj = reshape(uproj, [], 120, 3);
amps = log(1+sum(uproj.^2,3));
% amps = sum(uproj.^2,3)';
tall = ts;

%
[imax, amax] = max_interpolate(amps', rez, 50);
% hist(imax, 1000)

[tall, isort] = sort(tall);
amps = amps(isort);
imax = imax(isort);
amax = amax(isort);

mts = zeros(length(indBatch), 1);
for i = 1:length(indBatch)
   mts(i,1) = mean(tall(indBatch{i})); 
end
%% assign uproj spikes to clustered spikes
rezAM = rez;
rezAM.st3(:,2) = rez.st3(:,5);
rezAM.st3(rezAM.st3(:,2)==0, :) = []; 
[iclust, pmin] = match_uproj_to_rez(rezAM, imax, tall);

mean(pmin)

%%
[ampsort, isort] = sort(amax, 'descend');
bins = round(linspace(1, numel(amax), 100)); 

amp_bins = cell(numel(bins)-1, 1);
for i = 1:numel(bins)-1
   amp_bins{i, 1} = isort(bins(i):bins(i+1)-1); 
end

%
addpath(genpath('D:\CODE\MariusBox\Primitives'))
default_figure;
% close all
figure('Position', [1 1 12 6])
for j = 1:numel(amp_bins)
    gray_level = (j/numel(amp_bins)).^1.5; 
    plot(tall(amp_bins{j})/2.5e4, imax(amp_bins{j}), '.', 'MarkerSize', 3, 'Color', [1 1 1] * gray_level)
    hold on
end
cmap = colormap('jet');
cmap = cmap(ceil(64 * rand(max(iclust))), :);

if 1
    for j = 1:max(iclust)
            ix = (iclust==j) & (pmin<3);        
            plot(tall(ix)/2.5e4, imax(ix), '.', 'MarkerSize', 6, 'Color', cmap(j,:))
        hold on
    end
end



cmap = colormap('jet');
cmap = cmap(randperm(size(cmap,1)), :);
[uniqy] = unique(rez.yc);

hold on
for j = 1:numel(uniqy)
%    plot(mts/2.5e4, uniqy(j) + rez.totdY(j,:), 'Color', 'r', 'Linewidth', 1);
end

% xlim((0.5 + [.9 1.6])* 400)
% xlim([0 900]); ylim([700 1100]);
xlim([210 670]); ylim([860 1060]);

box off
ylabel('height (um)')
xlabel('Time (s)')
set(gcf, 'Color', 'w')

cd('D:\CODE\MariusBox\KilosortNIPS\fig')
%  print -dmeta drift_plot_with_lines.emf
export_fig('drift_plot_corrected_automerge_N192_zoomin.png')

%%