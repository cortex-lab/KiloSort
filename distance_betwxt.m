function [d2d, iY, drez] = distance_betwxt(W, U, mu, nbins)

Nfilt = numel(mu);
d2d = zeros(Nfilt);
for i = 1:size(U,3)
   d2d = d2d + (U(:,:,i)' * U(:,:,i)) .* (W(:,:,i)' * W(:,:,i));
end

mu = mu(:);
muall = repmat(mu, 1, Nfilt);
mrat = muall'./muall;
d2d = 1 - 2 * d2d./(1e-30 + mrat + 1./mrat);

d2d  = 1- triu(1 - d2d, 1);

% nbins = nbins(:)';
% nu = repmat(nbins, Nfilt, 1);

% d2d = d2d ./(1./nu + 1./nu');
% d2d = d2d ./(nu + nu');

[dMin, iY] = min(d2d, [], 1);

drez = dMin;

% drez = dMin .* nbins;

