
% iMrg = compute_merger(UW, cProj, st3, mu);
Cmerge = merge_posthoc(rez, mu, ops.fracse);
%%
vld = single(rez.nbins(1:Nfilt)>100);
links = (Cmerge < 5) & (vld * vld');

[nComponents,sizes,members] = networkComponents(links);

% nComponents

iMega = zeros(Nfilt, 1);
for i = 1:length(members)
   iMega(members{i}) = i; 
end
%

st3(:,5) = iMega(st3(:,2));


% members{iMega(126)}
nMega = zeros(nComponents, 1);
for j = 1:nComponents
    nMega(j) = sum(st3(:, 5)==j);
end
fprintf('Final number of merged clusters = %d \n', sum(nMega>100));
% %%
% figure;
% plot(U(:,[179 490 297],1))
% members{iMega(490)}
