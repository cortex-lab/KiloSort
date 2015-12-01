function [d2d, iY, drez] = distance_betwxt(W, U, mu, nbins)

Nfilt = numel(mu);
d2d = zeros(Nfilt);
for i = 1:size(U,3)
   d2d = d2d + (U(:,:,i)' * U(:,:,i)) .* (W(:,:,i)' * W(:,:,i));
end

mu = mu(:);
muall = repmat(mu, 1, Nfilt);
d2d = muall.^2 + muall'.^2 - 2 * muall.*muall'.*d2d;

d2d  = d2d + diag(Inf*ones(Nfilt,1));
[dMin, iY] = min(d2d, [], 1);

nbins = nbins(:)';
drez = dMin .* nbins;

