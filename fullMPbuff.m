Params(3) = 2;
% Params = double([NT Nfilt 0 10*maxFR 10]);

U0 = gpuArray(U);
utu = U0' * U0;
WtW = mexWtW(Params, W, utu);
WtW = permute(WtW, [3 1 2]);
rez.WtW = WtW;

clear nspikes
st3 = [];

if ops.verbose
   fprintf('Time %3.0fs. Running the final template matching pass...\n', toc) 
end

fid = fopen(fullfile(root, fnameTW), 'r');

for ibatch = 1:Nbatch
    if ibatch>Nbatch_buff
        offset = 2 * ops.Nchan*batchstart(ibatch-Nbatch_buff);
        fseek(fid, offset, 'bof');
        dat = fread(fid, [NT ops.Nchan], '*int16');
    else
       dat = DATA(:,:,ibatch); 
    end
    dataRAW = gpuArray(dat);
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    data 	= dataRAW * U; 
    
    [drez, dW, dU, st, id, x] = mexMPsub(Params,dataRAW,W,U,data,WtW);    
   
    nspikes(1:size(W,2)+1, ibatch) = histc(id, 0:1:size(W,2));
    delta(ibatch) = sum(x.^2)/1e6;
    STT = cat(2, double(st) +(NT-ops.ntbuff)*(ibatch-1), double(id)+1, double(x));
    st3 = cat(1, st3, STT);  
end
nsort = sort(sum(nspikes,2), 'descend');
fprintf('Time %3.0fs. ExpVar %2.6f, n10 %d, n20 %d, n30 %d, n40 %d \n', toc, nanmean(delta), nsort(10), nsort(20), ...
    nsort(min(size(W,2), 30)), nsort(min(size(W,2), 40)));

fclose(fid);
%%
fprintf('Time %3.0fs. Thresholding spikes at false positive rate...\n', toc) 
st3pos = [];
fprate = ops.fprate;
Thx = zeros(Nfilt,1);
for idd = 1:1:Nfilt
    ix = find(st3(:,2)==idd);
    xs = st3(ix, 3);
    
    Mu = 10*ops.Th;
    Nbins = 1000;
    
    bbins = linspace(0, Mu, Nbins);
    hpos = cumsum(hist(Mu - xs(xs>0), bbins));
    hneg = cumsum(hist(Mu + xs(xs<0), bbins));
    
    ifirst = find(hneg./hpos > fprate, 1);
    if isempty(ifirst)
        ifirst = numel(bbins);
    end
    Thx(idd) = Mu - bbins(ifirst);
    
    st3pos = cat(1, st3pos, st3(ix(xs>Thx(idd)), :));
end

[~, isort] = sort(st3pos(:,1), 'ascend');
st3pos = st3pos(isort,:);

rez.st3      = st3; 
rez.st3pos   = st3pos; 
rez.ops      = ops;

% WUnorms = sum(sum(dWUtotCPU.^2, 2), 1).^.5;
% rez.template = gather(dWUtotCPU ./ repmat(WUnorms, nt0, Nchan, 1));

rez.W = W;
rez.U = U;

for i = 1:Nfilt
    wav0 = W(:,i);
    wav0 = my_conv(wav0', .5)';
   [~, itrough] = min(wav0);
    [~, t2p] = max(wav0(itrough:end));
    rez.t2p(i,1) = t2p;
    rez.t2p(i,2) = itrough;   
end

rez.nbins = histc(rez.st3pos(:,2), .5:1:Nfilt+1);

[~, rez.ypos] = max(rez.U, [], 1);

%%
isV1pyr = rez.nbins(1:Nfilt)> 1000 & rez.ypos'>60 & rez.t2p(:,1)>10;
isV1pv = rez.nbins(1:Nfilt)> 1000 & rez.ypos'>60 & rez.t2p(:,1)<=10;
W0 = alignW(W);

figure(1)
hist(rez.t2p(rez.nbins>1000,1), 1:1:nt0) 
xlabel('trough to peak of templates with >1000 spikes')
set(gcf, 'Color', 'w')
ylabel('number of templates')
export_fig('fig1.pdf')

figure(2)
which_cells = find(isV1pyr);
[~, isort] = sort(rez.ypos(which_cells), 'ascend');
which_cells = which_cells(isort);
subplot(1,2,1)
imagesc(U(:, which_cells))
title('RS templates: spatial profile')
colormap('gray')
subplot(1,2,2)
plot(W0(:, which_cells))
axis tight
title('temporal profile')
set(gcf, 'Color', 'w')
export_fig('fig2.pdf')

figure(3)
which_cells = find(isV1pv);
[~, isort] = sort(rez.ypos(which_cells), 'ascend');
which_cells = which_cells(isort);
subplot(1,2,1)
imagesc(U(:, which_cells))
colormap('gray')

title('FS templates: spatial profile')
subplot(1,2,2)
plot(W0(:, which_cells))
axis tight
title('temporal profile')
set(gcf, 'Color', 'w')

export_fig('fig3.pdf')
%%

