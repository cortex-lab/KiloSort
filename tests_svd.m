igood   = logical(rez.st3pos(:,4));
coef    = rez.st3pos(igood, 5:end);
id      = rez.st3pos(igood, 2);

igood0   = logical(rez0.st3pos(:,4));
coef0    = rez0.st3pos(igood0, 5:end);
id0      = rez0.st3pos(igood0, 2);
%%
iNN = 187;
% iNN = iNN+1;
subplot(1,2,1)
plot(coef(id==iNN, 1), coef(id==iNN, 2), '.')
subplot(1,2,2)
plot(coef0(id0==iNN, 1), coef0(id0==iNN, 2), '.')

% plot(coef(id==iNN, 1), x(id==iNN), 'o')
%%
hist(rez.st3pos(igood,3), 1000)
%%
plot(sum(coefs,2), x, 'o')
%%

clear wu
for k = 1:Nfilt
    wu(:,:,k) = squeeze(W(:,k,:)) * squeeze(U(:,k,:))';

end