function [A,sWeight,B,Ekeep,Ob_keep,Ob_den_up,Ob_den_down] = do2TDVP_MPO(F,A,ML,M,MR,Band,OPTS,Site,OperatorR,OperatorL)
C_number = [0, 0; 0, 1];I = eye(2);


A0=A;
for i=1:Site-1
    A{i}=ncon({A{i},F},{[-1,1,-3],[-2,1]});
end
A{Site}=ncon({A{Site},OperatorR},{[-1,1,-3],[-2,1]});
%%%%% left-to-right 'warmup', put MPS in right orthogonal form
Nsites = length(A);
L{1} = ML; R{Nsites} = MR;
chid = size(M{1},3);
Mode_A=1;
for p = 1:Nsites - 1
    chil = size(A{p},1); chir = size(A{p},3);
    [qtemp,rtemp] = qr(reshape(A{p},[chil*chid,chir]),0);
    A{p} = reshape(qtemp,[chil,chid,chir]);
    Mode_A=Mode_A*norm(rtemp(:));
    A{p+1} = ncon({rtemp,A{p+1}},{[-1,1],[1,-2,-3]})/norm(rtemp(:));
    L{p+1} = ncon({L{p},M{p},A{p},conj(A{p})},{[2,1,4],[2,-1,3,5],[4,5,-3],[1,3,-2]});
end
chil = size(A{Nsites},1); chir = size(A{Nsites},3);
[qtemp,stemp] = qr(reshape(A{Nsites},[chil*chid,chir]),0);
A{Nsites} = reshape(qtemp,[chil,chid,chir]);
Mode_A=Mode_A*sqrt(trace(stemp*stemp'));% orthogonal constant
sWeight{Nsites+1} = stemp./sqrt(trace(stemp*stemp'));
W{Nsites} = reshape(qtemp*sWeight{Nsites+1},[chil,chid,chir]);

%initial state and spin flip

Ekeep =zeros(OPTS.numsweeps,1);
Ob_keep=zeros(OPTS.numsweeps,Nsites);
Ob_den_up=zeros(OPTS.numsweeps,Nsites);
Ob_den_down=zeros(OPTS.numsweeps,Nsites);
chimax=size(A{floor(Nsites/2)},1);
%%% Main program
for k = 1:OPTS.numsweeps
    for l=1:OPTS.midsweeps
        error=0;
        Mode_New=1;
        %%%%%% time_update
        %%%%%% Optimization sweep: right-to-left 
        for p = Nsites-1:-1:1
            %%%%% two-site update
            chil = size(A{p},1); chir = size(W{p+1},3);Length=chil*chid^2*chir;
            psi_t = reshape(ncon({A{p},W{p+1}},{[-1,-2,1],[1,-3,-4]}),[Length,1]);
            %%%%% Use Krylov method to update
            [psi_t] = Evol_Arnoldi(psi_t,min(OPTS.krydim,Length),-0.5j*OPTS.tau,@doApply_twosite,{L{p},M{p},M{p+1},R{p+1}});
%             [psi_t] = Evol_Lanczos(psi_t,min(OPTS.krydim,Length),-0.5j*OPTS.tau,@doApply_twosite,{L{p},M{p},M{p+1},R{p+1}});
            [utemp,stemp,vtemp] = svd(reshape(psi_t,[chil*chid,chid*chir]),'econ');
            chitemp = min(min(size(stemp)),chimax);
            A{p} = reshape(utemp(:,1:chitemp),[chil,chid,chitemp]);
            Mode=sqrt(sum(diag(stemp(1:chitemp,1:chitemp)).^2));Mode_all=sqrt(sum(diag(stemp).^2));
            error=error+abs(Mode_all-Mode)/Mode_all;
            Mode_New=Mode_New*Mode;
            sWeight{p+1} = stemp(1:chitemp,1:chitemp)./Mode;
            W{p}= ncon({A{p},sWeight{p+1}},{[-1,-2,1],[1,-3]});
            B{p+1} = reshape(vtemp(:,1:chitemp)',[chitemp,chid,chir]);

            %%%%% new block Hamiltonian MPO
            R{p} = ncon({M{p+1},R{p+1},B{p+1},conj(B{p+1})},{[-1,2,3,5],[2,1,4],[-3,5,4],[-2,3,1]});

            if p~=1
                chil = size(W{p},1); chir = size(W{p},3);Length=chil*chid*chir;
                psi_t = reshape(W{p},[Length,1]);
                %%%%% Use Krylov method to update
                [psi_t] = Evol_Arnoldi(psi_t,min(OPTS.krydim,Length),0.5j*OPTS.tau,@doApply_onesite,{L{p},M{p},R{p}});
%                 [psi_t] = Evol_Lanczos(psi_t,min(OPTS.krydim,Length),0.5j*OPTS.tau,@doApply_onesite,{L{p},M{p},R{p}});
                W{p} = reshape(psi_t,[chil,chid,chir]);
            end
        end

        %%%%%% left boundary tensor
        chil = size(W{1},1); chir = size(W{1},3);
        [utemp,stemp,vtemp] = svd(reshape(W{1},[chil,chid*chir]),'econ');
        B{1} = reshape(vtemp',[chil,chid,chir]);
        Mode=sqrt(trace(stemp.^2));
        Mode_New=Mode_New*Mode;
        sWeight{1} = utemp*stemp./Mode;

        %%%%%% Optimization sweep: left-to-right
        for p = 1:Nsites-1
            %%%%% two-site update
            chil = size(W{p},1); chir = size(B{p+1},3);Length=chil*chid^2*chir;
            psi_t = reshape(ncon({W{p},B{p+1}},{[-1,-2,1],[1,-3,-4]}),[Length,1]);
            %%%%% Use Krylov method to update
            [psi_t] = Evol_Arnoldi(psi_t,min(OPTS.krydim,Length),-0.5j*OPTS.tau,@doApply_twosite,{L{p},M{p},M{p+1},R{p+1}});
%             [psi_t] = Evol_Lanczos(psi_t,min(OPTS.krydim,Length),-0.5j*OPTS.tau,@doApply_twosite,{L{p},M{p},M{p+1},R{p+1}});
    
            [utemp,stemp,vtemp] = svd(reshape(psi_t,[chil*chid,chid*chir]),'econ');
            chitemp = min(min(size(stemp)),chimax);
            A{p} = reshape(utemp(:,1:chitemp),[chil,chid,chitemp]);
            Mode=sqrt(sum(diag(stemp(1:chitemp,1:chitemp)).^2));Mode_all=sqrt(sum(diag(stemp).^2));
            error=error+abs(Mode_all-Mode)/Mode_all;
            Mode_New=Mode_New*Mode;
            sWeight{p+1} = stemp(1:chitemp,1:chitemp)./Mode;
            B{p+1} = reshape(vtemp(:,1:chitemp)',[chitemp,chid,chir]);
            W{p+1}= ncon({sWeight{p+1},B{p+1}},{[-1,1],[1,-2,-3]});
                
            %%%%% new block Hamiltonian
            L{p+1} = ncon({L{p},M{p},A{p},conj(A{p})},{[2,1,4],[2,-1,3,5],[4,5,-3],[1,3,-2]});
    
            if p~=Nsites-1
                chil = size(W{p+1},1); chir = size(W{p+1},3);Length=chil*chid*chir;
                psi_t = reshape(W{p+1},[Length,1]);
                %%%%% Use Krylov method to update
                [psi_t] = Evol_Arnoldi(psi_t,min(OPTS.krydim,Length),0.5j*OPTS.tau,@doApply_onesite,{L{p+1},M{p+1},R{p+1}});
%                 [psi_t] = Evol_Lanczos(psi_t,min(OPTS.krydim,Length),0.5j*OPTS.tau,@doApply_onesite,{L{p+1},M{p+1},R{p+1}});
                W{p+1} = reshape(psi_t,[chil,chid,chir]);
            end
        end
        
        %%%%%% right boundary tensor
        chil = size(W{Nsites},1); chir = size(W{Nsites},3);
        [utemp,stemp,vtemp] = svd(reshape(W{Nsites},[chil*chid,chir,1]),'econ');    
        A{Nsites} = reshape(utemp,[chil,chid,chir]);
        Mode=sqrt(trace(stemp.^2));
        Mode_New=Mode_New*Mode;
        sWeight{Nsites+1} = (stemp./Mode)*vtemp';
        if error>1e-9
            chimax=min(chimax+Band.chi_step,Band.chimax);
        elseif error<1e-9
            chimax=max(chimax-Band.chi_step,Band.chimin);
        end
    end
    LL = ncon({L{Nsites},M{Nsites},W{Nsites},conj(W{Nsites})},{[2,1,4],[2,-1,3,5],[4,5,-3],[1,3,-2]});
    Ekeep(k,1)=abs(Mode_A)^2*ncon({LL,R{Nsites}},{[1,2,3],[1,2,3]});
    Ob_keep(k,:) = Calculate_innerP(F,A,A0,Mode_A,Nsites,OperatorL);
    Ob_den_up(k,:) = DMRG_OneSiteObservation(A,sWeight,kron(C_number,I));
    Ob_den_down(k,:) = DMRG_OneSiteObservation(A,sWeight,kron(I,C_number));
    fprintf('Sweep: %2.1d of %2.1d, CutEr:%2.2d, Mode:%6.6d, ReN:%6.6d, Energy: %12.12d, Band: %2.2d\n',k,OPTS.numsweeps,error,Mode_New,...
        sum(Ob_den_up(k,:)+Ob_den_down(k,:)),real(Ekeep(k,1)),chimax);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Result = Calculate_innerP(F,A,A0,sWeight,Nsites,OperatorL)
% applies the superblock MPO to the state
Result=zeros(1,Nsites);chil=size(A{1},1);chir=size(A{Nsites},3);
L_Psi=cell(1,Nsites+1);R_Psi=cell(1,Nsites+1);
L_Psi{1}=ones(chil,chil);
R_Psi{Nsites+1}=ones(chir,chir);
for p=1:Nsites
    L_Psi{p+1}=ncon({L_Psi{p},A{p},F,conj(A0{p})},{[1,4],[1,2,-1],[3,2],[4,3,-2]});
    R_Psi{Nsites-p+1}=ncon({A{Nsites-p+1},conj(A0{Nsites-p+1}),R_Psi{Nsites-p+2}},{[-1,2,1],[-2,2,3],[1,3]});
end

for p=1:Nsites
    Result(1,p)=sWeight*ncon({L_Psi{p},A{p},OperatorL,conj(A0{p}),R_Psi{p+1}},{[1,4],[1,2,5],[3,2],[4,3,6],[5,6]});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = doApply_twosite(psi,L,M1,M2,R)
% applies the superblock MPO to the state

psi = reshape(ncon({reshape(psi,[size(L,3),size(M1,4),size(M2,4),size(R,3)]),L,M1,M2,R},...
    {[1,3,5,7],[2,-1,1],[2,4,-2,3],[4,6,-3,5],[6,-4,7]}),[size(L,3)*size(M1,4)*size(M2,4)*size(R,3),1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = doApply_onesite(psi,L,M,R)
% applies the superblock MPO to the state
psi = reshape(ncon({reshape(psi,[size(L,3),size(M,4),size(R,3)]),L,M,R},...
    {[1,3,5],[2,-1,1],[2,4,-2,3],[4,-3,5]}),[size(L,2)*size(M,4)*size(R,2),1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psivec] = Evol_Arnoldi(psivec,krydim,dt,linFunct,functArgs)
% function for computing the smallest algebraic eigenvalue and eigenvector
% of the linear function 'linFunct' using a Lanczos method. Maximum
% iterations are specified by 'OPTS.maxit' and the dimension of Krylov
% space is specified by 'OPTS.krydim'. Input 'functArgs' is an array of
% 114514
% optional arguments passed to 'linFunct'.
if norm(psivec) == 0
    psivec = rand(length(psivec),1);
end
psi = zeros(numel(psivec),krydim+1);
A = zeros(krydim+1,krydim);
%Calculate krylov subspace
psi(:,1) = psivec(:)/norm(psivec);
for p = 2:krydim+1
    psi0 = linFunct(psi(:,p-1),functArgs{(1:length(functArgs))});
    psi(:,p) = psi0;
    for g = 1:1:p-1
        A(g,p-1)=dt*dot(psi(:,g),psi0);
    end
    for g = 1:1:p-1
        h=dot(psi(:,g),psi(:,p));
        psi(:,p) = psi(:,p) - h*psi(:,g);
        psi(:,p) = psi(:,p)/max(norm(psi(:,p)),1e-16);
    end
    A(p,p-1)=dt*dot(psi(:,p),psi0);
end

U_t=expm(A(1:krydim,1:krydim));
Psi_t=U_t(:,1)*norm(psivec);
psivec = psi(:,1:krydim)*Psi_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psivec] = Evol_Lanczos(psivec,krydim,dt,linFunct,functArgs)
% function for computing the smallest algebraic eigenvalue and eigenvector
% of the linear function 'linFunct' using a Lanczos method. Maximum
% iterations are specified by 'OPTS.maxit' and the dimension of Krylov
% space is specified by 'OPTS.krydim'. Input 'functArgs' is an array of
% optional arguments passed to 'linFunct'.

if norm(psivec) == 0
    psivec = rand(length(psivec),1);
end
psi = zeros(numel(psivec),krydim+1);
A = zeros(krydim,krydim);

psi(:,1) = psivec(:)/norm(psivec);
for p = 2:krydim+1
    psi(:,p) = linFunct(psi(:,p-1),functArgs{(1:length(functArgs))});
    for g = 1:1:p-1
        A(p-1,g) = dt*dot(psi(:,p),psi(:,g));
        A(g,p-1) = dt*conj(A(p-1,g));
    end
    for g = 1:1:p-1
        psi(:,p) = psi(:,p) - dot(psi(:,g),psi(:,p))*psi(:,g);
        psi(:,p) = psi(:,p)/max(norm(psi(:,p)),1e-16);
    end
end

U_t=expm(A(1:krydim,1:krydim));
Psi_t=U_t(:,1)*norm(psivec);
psivec = psi(:,1:krydim)*Psi_t;