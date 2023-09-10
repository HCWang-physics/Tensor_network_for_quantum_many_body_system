OPTS.maxit = 100; % iterations of Lanczos method
OPTS.krydim = 6; % dimension of Krylov subspace

L=1000;t=1;gamma=0;
psivec=rand(L,1);
H=diag((t-gamma)*ones(1,L-1),1)+diag((t+gamma)*ones(1,L-1),-1)+diag(0j*ones(1,L),0);

[psivec,dval,A] = eig_Lancos(psivec,H,OPTS);
[Value]=eig(H);
[Y,I]=sort(real(Value));
Value(I(1))
%%
figure(1)
subplot(121)
plot(real(dval))
subplot(122)
plot(imag(dval))
