function out = get_dySpikes(Mmax)

Mmax = bsxfun(@minus, Mmax, max(Mmax, [], 2));

ix = (Mmax(:,1) > Mmax(:,2)) & (Mmax(:,3) > Mmax(:,2));

out = (Mmax(:,3) - Mmax(:,1))./ (2 * max(abs(Mmax(:,1)), abs(Mmax(:,3))));

out(ix) = 0;
out = min(1, max(-1, out));