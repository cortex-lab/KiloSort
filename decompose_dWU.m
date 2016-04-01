function  [W, U, mu, UtU] = decompose_dWU(dWUtot, Nrank)

[nt0 Nchan Nfilt] = size(dWUtot);

W = zeros(nt0, Nrank, Nfilt, 'single');
U = zeros(Nchan, Nrank, Nfilt, 'single');
mu = zeros(Nfilt, 1, 'single');

dWUtot(isnan(dWUtot)) = 0;
for k = 1:Nfilt
    
    [Wall, Sv, Uall] = svd(gather(dWUtot(:,:,k)), 0);
    Wall = Wall * Sv; 
    
    Sv = diag(Sv);
    mu(k) = sum(Sv(1:Nrank).^2).^.5;
    Wall = Wall/mu(k);
    [~, imax] = max(abs(Wall(:,1)));
    Wall(:,1) = -Wall(:,1) * sign(Wall(imax,1));
    
    W(:,:,k) = Wall(:,1:Nrank);
    U(:,:,k) = Uall(:,1:Nrank);

end
U = permute(U, [1 3 2]);
W = permute(W, [1 3 2]);

U(isnan(U)) = 0;

UtU = abs(U(:,:,1)' * U(:,:,1)) > .1;

% mu = min(mu, 200);