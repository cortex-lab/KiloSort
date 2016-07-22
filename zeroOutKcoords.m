function U = zeroOutKcoords(U, kcoords)

[M, imax] = squeeze(max(abs(U(:,:,1)), [], 1));

for i = 1:size(U,2)
   U(kcoords~=kcoords(imax(i)),i,:) = 0;
   U(:,i,:) = normc(squeeze(U(:,i,:)));
end

