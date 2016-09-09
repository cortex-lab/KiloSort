function [alldy, allNSpikes] = bin_corrections(dySpikes, ySp, Nmax)

alldy = zeros(Nmax, 1);
allNSpikes = zeros(Nmax, 1);
for i = 1:Nmax
    ifind = find(ySp==i);
    if ~isempty(ifind)
        alldy(i) = sum(dySpikes(ifind));
        allNSpikes(i) = numel(ifind);
    end
end

