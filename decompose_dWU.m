function  [W, U, mu, UtU, nu] = decompose_dWU(dWU, Nrank)

[nt0 Nchan Nfilt] = size(dWU);

W = zeros(nt0, Nrank, Nfilt, 'single');
U = zeros(Nchan, Nrank, Nfilt, 'single');
mu = zeros(Nfilt, 1, 'single');
% dmax = zeros(Nfilt, 1);

dWU(isnan(dWU)) = 0;
for k = 1:Nfilt
    
    [Wall, Sv, Uall] = svd(gather(dWU(:,:,k)), 0);
    [~, imax] = max(abs(Wall(:,1)));
    Uall(:,1) = -Uall(:,1) * sign(Wall(imax,1));
    Wall(:,1) = -Wall(:,1) * sign(Wall(imax,1));
     
%     [~, imin] = min(diff(Wall(:,1), 1));
%     [~, imin] = min(Wall(:,1));
%     dmax(k) = - (imin- 20);
    
%     if dmax(k)>0
%         dWU((dmax(k) + 1):nt0, :,k) = dWU(1:nt0-dmax(k),:, k);
%         Wall((dmax(k) + 1):nt0, :)  = Wall(1:nt0-dmax(k),:);
%     else
%         dWU(1:nt0+dmax(k),:, k) = dWU((1-dmax(k)):nt0,:, k);
%         Wall(1:nt0+dmax(k),:) = Wall((1-dmax(k)):nt0,:);
%     end
    
    Wall = Wall * Sv; 
    
    Sv = diag(Sv);
    mu(k) = sum(Sv(1:Nrank).^2).^.5;
    Wall = Wall/mu(k);
   
    W(:,:,k) = Wall(:,1:Nrank);
    U(:,:,k) = Uall(:,1:Nrank);

end
U = permute(U, [1 3 2]);
W = permute(W, [1 3 2]);

U(isnan(U)) = 0;

UtU = abs(U(:,:,1)' * U(:,:,1)) > .1;


Wdiff = cat(1, W, zeros(2, Nfilt, Nrank)) - cat(1,  zeros(2, Nfilt, Nrank), W);
nu = sum(sum(Wdiff.^2,1),3);
nu = nu(:);



% mu = min(mu, 200);