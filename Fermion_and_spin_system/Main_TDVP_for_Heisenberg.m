% mainDMRG_MPO
% ------------------------ 
tic

%%%%% Example 1: Heisenberg model %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Set simulation options
chi = 100; % maximum bond dimension
Nsites = 100; % number of lattice sites

%%%% Define Hamiltonian MPO (quantum XX model)
chid = 2;J=0.3;Delta=1/J;hm=0.1; 
Sx = 1/2 * [0, 1; 1, 0];
Sy = 1/2 * [0, -1i; 1i, 0];
Sz = 1/2 * [1, 0; 0, -1];
sI = eye(2);
Splus = Sx + 1i* Sy;
H0=hm*Sz;

M0 = zeros(5,5,2,2);
M0(1,1,:,:) = sI; M0(1,5,:,:) = H0;
M0(1,2,:,:) = sqrt(0.5*J)*Splus; M0(2,5,:,:) = sqrt(0.5*J)*Splus';
M0(1,3,:,:) = sqrt(J*Delta)*Sz; M0(3,5,:,:) = sqrt(J*Delta)*Sz;
M0(1,4,:,:) = sqrt(0.5*J)*Splus'; M0(4,5,:,:) = sqrt(0.5*J)*Splus;
M0(5,5,:,:) = sI;
ML = reshape([1;0;0;0;0],[5,1,1]); %left MPO boundary
MR = reshape([0;0;0;0;1],[5,1,1]); %right MPO boundary

%%%% Initialize MPS tensors
A_initial = {};M = {};
A_initial{1} = rand(1,chid,min(chi,chid));
M{1}=M0;
for k = 2:Nsites
    A_initial{k} = rand(size(A_initial{k-1},3),chid,min(min(chi,size(A_initial{k-1},3)*chid),chid^(Nsites-k)));
    M0(1,5,:,:)=-M0(1,5,:,:);
    M{k}=M0;
end

%%
%%%% Do DMRG sweeps
OPTS.numsweeps = 4; % number of DMRG sweeps
OPTS.display = 1; % level of output display
OPTS.updateon = 1; % update methond 1=Arnoldi 2=eigLanczos
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 2; % dimension of Krylov subspace
[A0,sWeight0,B0,Ekeep0] = doDMRG_MPO(A_initial,ML,M,MR,chi,OPTS);
%%
TDVP.numsweeps = 4000; % number of time iteration
TDVP.midsweeps = 2; % number of time iteration
TDVP.tau = 0.01; % time step
TDVP.krydim=6; % dimension of Krylov subspace
[A,sWeight,B,Ekeep1,Ob_keep] = do2TDVP_MPO(A0,ML,M,MR,chi,TDVP,floor(Nsites/2),Sz);
%%
[Twosite_correlation] = DMRG_TwoSiteObservation(A0,sWeight0,Sz,floor(Nsites/2));
Szz=zeros(TDVP.numsweeps,Nsites);
for i=1:TDVP.numsweeps
    T=TDVP.tau*(i)*TDVP.midsweeps;
    Szz(i,:)=exp(1j*Ekeep0(end)*T)*Ob_keep(i,:)-Twosite_correlation;
end
%%
[Ob_M] = DMRG_OneSiteObservation(A0,sWeight0,Sz);
[Twosite_correlation] = DMRG_TwoSiteObservation(A0,sWeight0,Sz,floor(Nsites/2));
%%
figure(1);
subplot(121)
plot(1:Nsites,real(Ob_M), '-o', 'disp', 'Tensor network DMRG')
xlabel('i')
ylabel('<S_{z}>')
subplot(122)
plot(1:Nsites,real(Twosite_correlation), '-mh')
xlabel('|i-j|')
ylabel('correlation')
%%
Time=1:TDVP.numsweeps;
Time=Time*TDVP.tau*TDVP.midsweeps;
[X,Y]=meshgrid(1:Nsites,Time);

Number=100; 
Fourier0=zeros(Number,Nsites);
Fourier=zeros(Number,Nsites+1);
Omega=1:Number;
Omega=Omega*2*pi/(TDVP.numsweeps*TDVP.midsweeps*TDVP.tau);
for i=1:Number
    Time_exp=exp(-1j*Time*Omega(i));
    Fourier0(i,:)=sum(Time_exp'.*real(Szz),1)/(TDVP.numsweeps*TDVP.midsweeps*TDVP.tau);
end

KK_x=0:Nsites;
KK_x=KK_x*2*pi/Nsites;
XX=1:Nsites;
for i=1:Nsites+1
    Real_exp=exp(1j*XX*KK_x(i));
    Fourier(:,i)=sum(Real_exp.*Fourier0,2)/Nsites;
end
%%
[X0,Y0]=meshgrid(KK_x,Omega);
figure(2);
subplot(121)
pcolor(X0/pi,Y0/J,abs(Fourier));
colorbar
colormap('hot')
% caxis([-0.5,0.5])
shading interp
xlim([0,2])
xlabel('k_{x}(\pi)')
ylabel('E(J)')
subplot(122)
pcolor(X,Y*J,log(abs(Szz)));
colorbar
caxis([-10,-2])
shading interp
xlabel('Site')
ylabel('T(1/J)')

%%
toc