function Cmerge = merge_posthoc(rez, mu, fracse)
%fracse = 0.1;
ops = rez.ops;
LAM = ops.lam(3) * (20./mu).^2;
Nfilt = rez.ops.Nfilt;

Cmerge = Inf *ones(Nfilt);

for testID = 1:Nfilt
    clusterIDs = rez.st3(:,2);
    tfi = rez.iNeigh;
    tf = rez.cProj;
    
    spikesTest = find(clusterIDs==testID);
    
    simIDs = tfi(:,testID);
        
    if numel(spikesTest)>20 
        for s = 1:length(simIDs)
            simS_T = find(tfi(:,simIDs(s))==testID);
            spikesSim = find(clusterIDs==simIDs(s));
            
            if simIDs(s)>testID && numel(spikesSim)>200 && ~isempty(simS_T)
                ft1 = [tf(spikesTest,1); tf(spikesSim,simS_T)];
                ft2 = [tf(spikesTest,s); tf(spikesSim,1)];
                
                ft1 = (ft1 + LAM(testID) * mu(testID))    / sqrt(1 + LAM(testID));
                ft2 = (ft2 + LAM(simIDs(s)) * mu(simIDs(s))) / sqrt(1 +  LAM(simIDs(s)));
                
                df = ft1 - ft2;
                l1 = min(df(:));
                l2 = max(df(:));
                
                df1 = df(1:numel(spikesTest));
                df2 = df(1+numel(spikesTest):end);
                
                se = (std(df1) + std(df2))/2;
                se25 = fracse * se;
                b2 = [0:se25:-l1];
                b1 = [0:se25:l2];
                
                hs1 = my_conv(histc(df1, b1), 1);
                hs2 = my_conv(histc(-df2, b2), 1);
                
                mmax = min(max(hs1), max(hs2));
                
                m1 = ceil(mean(df1)/se25);
                m2 = -ceil(mean(df2)/se25);
                steps = sum(hs1(1:m1)<mmax/5) + sum(hs2(1:m2)<mmax/5);
                
%                 steps = find(hs1>mmax/5, 1) + find(hs2>mmax/5, 1) - 2;
                Cmerge(testID, simIDs(s)) =steps;
                
                %                 mlow = min(max(hs1(1), hs1(2)), max(hs2(1), hs2(2)));
                %                 Cmerge(testID, simIDs(s)) = min( max(hs2)/mlow, max(hs1)/mlow);
                Cmerge(simIDs(s),testID) = Cmerge(testID, simIDs(s));
            end
        end
    end
end
