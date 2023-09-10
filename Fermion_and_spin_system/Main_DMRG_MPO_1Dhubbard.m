% mainDMRG_MPO
tic
%%%%% Set simulation options
chi = 200; % maximum bond dimension
Nsites = 20; % number of lattice sites

%%%% Define Hamiltonian MPO (quantum XX model)
chid = 4;J=exp(-0);Delta=0;U=4;hm=0;
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
A_initial{1} = rand(1,chid,min(chi,chid));
M_initial{1}=M0;
for k = 2:Nsites
    A_initial{k} = rand(size(A_initial{k-1},3),chid,min(min(chi,size(A_initial{k-1},3)*chid),chid^(Nsites-k)));
    M_initial{k}=M0;
end

M=M_initial;
%[M] = Normalize_MPO(M_initial);
%%
%%%% Do DMRG sweeps
OPTS.numsweeps = 10; % number of DMRG sweeps
OPTS.display = 1; % level of output display
OPTS.updateon = 1; % update methond 1=Arnoldi 2=eigLanczos
OPTS.maxit = 2; % iterations of Lanczos method
OPTS.krydim = 2; % dimension of Krylov subspace
[A,sWeight,B,Ekeep] = doDMRG_MPO(A_initial,ML,M,MR,chi,OPTS);
%%
[Ob_M_up] = DMRG_OneSiteObservation(A,sWeight,kron(C_number,I));
[Ob_M_down] = DMRG_OneSiteObservation(A,sWeight,kron(I,C_number));
[Ob_Sz] = DMRG_OneSiteObservation(A,sWeight,Sz);
[Twosite_correlation_up] = DMRG_TwoSiteObservation(F,A,sWeight,up_plus',up_plus,1,floor(1));
[Twosite_correlation_down] = DMRG_TwoSiteObservation(F,A,sWeight,down_plus',down_plus,1,floor(1));
[Twosite_correlation_Sz] = DMRG_TwoSiteObservation(F,A,sWeight,Sz,Sz,2,floor(1));
%%
figure(1);
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
ylabel('<n_{down}>')
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
toc