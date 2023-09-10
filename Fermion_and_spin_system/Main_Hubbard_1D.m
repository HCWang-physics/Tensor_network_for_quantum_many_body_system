% mainDMRG_MPO
% ------------------------ 
clear
tic
%%%%% Example 1: 1D hurbbard model (Green function) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Set simulation options
Band.chimin = 2^6; % maximum bond dimension
Band.chi_step=2^6;
Band.chimax=2^8;
Nsites = 20; % number of lattice sites

%%%% Define Hamiltonian MPO (quantum XX model)
chid = 4;J=-0.5;Delta=0;U=1;hm=0;
C_plus = [0, 0; 1, 0];
C_number = [0, 0; 0, 1];
I = eye(2);
F=[1,0,0,0;0,-1,0,0;0,0,-1,0;0,0,0,1];
up_plus=kron(C_plus,I);
down_plus=kron(I,C_plus);
down_plus=F*down_plus;
SI=eye(4);Sz=0.5*(kron(C_number,I)-kron(I,C_number));
Sx=0.5*(kron(C_plus,C_plus')+kron(C_plus',C_plus));
Sy=0.5j*(-1*kron(C_plus,C_plus')+kron(C_plus',C_plus));
H0=U*kron(C_number-0.5*eye(2),C_number-0.5*eye(2))+Delta*kron(C_number,I)+Delta*kron(I,C_number)+hm*Sz;

M0 = zeros(6,6,chid,chid);
M0(1,1,:,:) = SI; M0(1,6,:,:) = H0;
M0(1,2,:,:) = J*up_plus*F; M0(2,6,:,:) = up_plus';
M0(1,3,:,:) = J*down_plus*F; M0(3,6,:,:) = down_plus';
M0(1,4,:,:) = -J*down_plus'*F; M0(4,6,:,:) = down_plus;
M0(1,5,:,:) = -J*up_plus'*F; M0(5,6,:,:) = up_plus;
M0(6,6,:,:) = SI;
ML = reshape([1;0;0;0;0;0],[6,1,1]); %left MPO boundary
MR = reshape([0;0;0;0;0;1],[6,1,1]); %right MPO boundary

%%%% Initialize MPS tensors
A_initial = {};M_initial = {};
A_initial{1} = rand(1,chid,min(Band.chimin,chid));
M_initial{1}=M0;
for k = 2:Nsites
    A_initial{k} = rand(size(A_initial{k-1},3),chid,min(min(Band.chimin,size(A_initial{k-1},3)*chid),chid^(Nsites-k)));
    M_initial{k}=M0;
end

M=M_initial;
%[M] = Normalize_MPO(M_initial);
%%
%%%% Do DMRG sweeps
OPTS.numsweeps = 20; % number of DMRG sweeps
OPTS.display = 1; % level of output display
OPTS.updateon = 1; % update methond 1=Arnoldi 2=eigLanczos
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 4; % dimension of Krylov subspace
[A0,sWeight0,B0,Ekeep0] = doDMRG_MPO(A_initial,ML,M,MR,Band,OPTS);
%%
% [Ob_M_up0] = DMRG_OneSiteObservation(A0,sWeight0,kron(C_number,I));
% [Ob_M_down0] = DMRG_OneSiteObservation(A0,sWeight0,kron(I,C_number));
% [Ob_Sz0] = DMRG_OneSiteObservation(A0,sWeight0,Sz);
% [Twosite_correlation_up0] = DMRG_TwoSiteObservation(F,A0,sWeight0,up_plus',up_plus,1,floor(Nsites/2));
% [Twosite_correlation_down0] = DMRG_TwoSiteObservation(F,A0,sWeight0,down_plus',down_plus,1,floor(Nsites/2));
% [Twosite_correlation_Sz0] = DMRG_TwoSiteObservation(F,A0,sWeight0,Sz,Sz,2,floor(Nsites/2));

% figure(1);
% subplot(411)
% Site_L=length(Ekeep0);
% plot(1:Site_L,real(Ekeep0)/Nsites,'-k')
% hold on
% plot(1:Site_L,imag(Ekeep0)/Nsites,'--r')
% xlabel('Update step')
% ylabel('Ground energy')
% legend('Re{E0}','Im{E0}')
% subplot(423)
% plot(1:Nsites,real(Ob_M_up0), '-o', 'disp', 'Tensor network DMRG')
% xlabel('i')
% ylabel('<n_{up}>')
% subplot(424)
% plot(1:Nsites,real(Ob_M_down0), '-o', 'disp', 'Tensor network DMRG')
% xlabel('i')
% ylabel('<n_{down}>')
% subplot(425)
% plot(1:Nsites,real(Ob_M_up0)+real(Ob_M_down0), '-o', 'disp', 'Tensor network DMRG')
% xlabel('i')
% ylabel('<n_{down}+n_{up}>')
% subplot(426)
% plot(1:Nsites,real(Ob_Sz0), '-o', 'disp', 'Tensor network DMRG')
% xlabel('i')
% ylabel('<S_{z}>')
% subplot(427)
% plot(1:Nsites,real(Twosite_correlation_up0), '-mh')
% hold on
% plot(1:Nsites,real(Twosite_correlation_down0), '-sr')
% xlabel('|i-j|')
% legend('Cup^{+}Cup','Cdown^{+}Cdown')
% ylabel('correlation')
% ylim([-0.25,0.5])
% subplot(428)
% plot(1:Nsites,real(Twosite_correlation_Sz0), '-ob')
% xlabel('|i-j|')
% legend('SzSz')
% ylabel('correlation')
% ylim([-0.2,0.25])
%%
TDVP.numsweeps = 4000; % number of time iteration
TDVP.midsweeps = 2; % number of time iteration
TDVP.tau = 0.005; % time step
TDVP.krydim=4; % dimension of Krylov subspace
[A,sWeight,B,Ekeep,Ob_keep,Ob_den_up,Ob_den_down] = do2TDVP_MPO(F,A0,ML,M,MR,Band,TDVP,floor(Nsites/2),Sz,Sz);
%%
[Twosite_correlation] = DMRG_TwoSiteObservation(F,A0,sWeight0,Sz,Sz,2,floor(Nsites/2));
Coree=zeros(TDVP.numsweeps,Nsites);
for i=1:TDVP.numsweeps
    T=TDVP.tau*(i)*TDVP.midsweeps;
    Coree(i,:)=exp(1j*Ekeep0(end)*T)*Ob_keep(i,:)-Twosite_correlation;
end
%%
[Ob_M_up] = DMRG_OneSiteObservation(A,sWeight,kron(C_number,I));
[Ob_M_down] = DMRG_OneSiteObservation(A,sWeight,kron(I,C_number));
[Ob_Sz] = DMRG_OneSiteObservation(A,sWeight,Sz);
[Twosite_correlation_up] = DMRG_TwoSiteObservation(F,A,sWeight,up_plus',up_plus,1,floor(1));
[Twosite_correlation_down] = DMRG_TwoSiteObservation(F,A,sWeight,down_plus',down_plus,1,floor(1));
[Twosite_correlation_Sz] = DMRG_TwoSiteObservation(F,A,sWeight,Sz,Sz,2,floor(1));
%%
figure(2);
subplot(411)
Site_L=length(Ekeep);
plot(1:Site_L,real(Ekeep)/Nsites,'-k')
hold on
plot(1:Site_L,imag(Ekeep)/Nsites,'--r')
xlabel('Update step')
ylabel('Ground energy')
legend('Re{E0}','Im{E0}')
subplot(423)
plot(1:Nsites,real(Ob_M_up), '-o', 'disp', 'Tensor network DMRG')
xlabel('i')
ylabel('<n_{up}>')
subplot(424)
plot(1:Nsites,real(Ob_M_down), '-o', 'disp', 'Tensor network DMRG')
xlabel('i')
ylabel('<n_{down}>')
subplot(425)
plot(1:Nsites,real(Ob_M_up)+real(Ob_M_down), '-o', 'disp', 'Tensor network DMRG')
xlabel('i')
ylabel('<n_{down}+n_{up}>')
subplot(426)
plot(1:Nsites,real(Ob_Sz), '-o', 'disp', 'Tensor network DMRG')
xlabel('i')
ylabel('<S_{z}>')
subplot(427)
plot(1:Nsites,real(Twosite_correlation_up), '-mh')
hold on
plot(1:Nsites,real(Twosite_correlation_down), '-sr')
xlabel('|i-j|')
legend('Cup^{+}Cup','Cdown^{+}Cdown')
ylabel('correlation')
ylim([-0.25,0.5])
subplot(428)
plot(1:Nsites,real(Twosite_correlation_Sz), '-ob')
xlabel('|i-j|')
legend('SzSz')
ylabel('correlation')
ylim([-0.2,0.25])
%%
Time=1:TDVP.numsweeps;
Time=Time*TDVP.tau*TDVP.midsweeps;
[X,Y]=meshgrid(1:Nsites,Time);

Number=60; 
Fourier0=zeros(Number,Nsites);
Fourier=zeros(Number,Nsites+1);
Omega=1:Number;
Omega=Omega*2*pi/(TDVP.numsweeps*TDVP.midsweeps*TDVP.tau);
for i=1:Number
    Time_exp=exp(-1j*Time*Omega(i));
    Fourier0(i,:)=sum(Time_exp'.*Coree,1)/(TDVP.numsweeps*TDVP.midsweeps*TDVP.tau);
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
figure(3);
subplot(131)
pcolor(X0/pi,Y0,abs(Fourier));
colormap('hot')
shading interp
xlabel('k_{x}(\pi)')
ylabel('E(U)')
subplot(132)
pcolor(X,Y,abs(Coree));
colorbar
shading interp
xlabel('Site')
ylabel('T(1/U)')
subplot(133)
pcolor(X,Y,abs(Ob_keep));
colorbar
shading interp
xlabel('Site')
ylabel('T(1/U)')

%%
toc