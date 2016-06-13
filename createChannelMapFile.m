%  create a channel map file

Nchannels = 32;
connected = ones(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

save('C:\DATA\Spikes\20150601_chan32_4_900s\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind')

%%

Nchannels = 32;
connected = ones(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords   = xcoords(:);
ycoords   = repmat(1:Nchannels/4, 4, 1);
ycoords   = ycoords(:);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

save('C:\DATA\Spikes\Piroska\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind')
