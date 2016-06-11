function [dWU, st, id, x,Cost, nsp] = ...
            mexMPregMUcpu(Params,dataRAW,fW,data,UtU,mu, lam , dWU, nu, ops);

NT      = Params(1);
nFilt   = Params(2);

pm      = Params(8);
fdata   = fft(data, [], 1);
proj       = real(ifft(fdata .* fW(:,:), [], 1));
proj       = sum(reshape(proj, NT, nFilt, 3),3);
%
Ci = bsxfun(@plus, proj, (mu.*lam)');
Ci = bsxfun(@rdivide, Ci.^2,  1 + lam');
Ci = bsxfun(@minus, Ci, (lam .*mu.^2)');
%
[mX, id] = max(Ci,[], 2);
maX     = -my_min(-mX, 21, 1);

st      = find(maX < mX + 1e-6);
st(st>NT-61) = [];

id      = id(st);

inds = bsxfun(@plus, st', [1:61]'-1);
dspk = reshape(dataRAW(inds, :), 61, numel(st), ops.Nchan);
dspk = permute(dspk, [1 3 2]);

x = zeros(size(id));
Cost = zeros(size(id));
nsp = zeros(nFilt,1);
for j = 1:size(dspk,2)
    dWU(:,:,id(j)) = pm * dWU(:,:,id(j)) + (1-pm) * dspk(:,:,j);
    x(j) = proj(st(j), id(j));
    Cost(j) = maX(st(j));
    nsp(id(j)) = nsp(id(j)) + 1;
end

id = id - 1;
