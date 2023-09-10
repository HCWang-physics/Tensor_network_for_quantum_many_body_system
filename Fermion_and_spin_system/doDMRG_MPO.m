function [A,sWeight,B,Ekeep] = doDMRG_MPO(A,ML,M,MR,Band,OPTS)
% function [A,sWeight,B,Ekeep] = doDMRG_MPO(A,ML,M,MR,chimax,OPTS)
% ------------------------ 
% by Glen Evenbly (c) for www.tensors.net, (v1.1) - last modified 21/1/2019
% 
% Implementation of DMRG for a 1D chain with open boundaries, using the
% two-site update strategy. Each update is accomplished using a custom
% implementation of the Lanczos iteration to find (an approximation to) the
% ground state of the superblock Hamiltonian. Input 'A' is containing the
% MPS tensors whose length is equal to that of the 1D lattice. The
% Hamiltonian is specified by an MPO with 'ML' and 'MR' the tensors at the
% left and right boundaries, and 'M' the bulk MPO tensor. Automatically
% grow the MPS bond dimension to maximum dimension 'chimax'. Outputs 'A'
% and 'B' are cells of the MPS tensors in left and right orthogonal form
% respectively, while 'sWeight' is a cell of the Schmidt coefficients
% across different lattice positions. 'Ekeep' is a vector describing the
% energy at each update step.
%
% Options:
% OPTS.numsweeps: number of DMRG sweeps [integer | {10}]
% OPTS.display: 0 => no display, 
%               1 => display at end of each sweep, 
%               2 => data displayed at each step, [0 | {1} | 2] 
% OPTS.updateon: 0 => MPS not updated (but is brought into cannonical form), 
%                1 => DMRG using 2-site updates
% OPTS.maxit: iterations of Lanczos method for each superblock diagonalization [integer | {2}]
% OPTS.krydim: dimension of Krylov space in superblock diagonalization [integer > 1, {4}]

%%%%% set options to defaults if not specified 
if OPTS.krydim < 2
    OPTS.krydim = 2;
end

%%%%% left-to-right 'warmup', put MPS in right orthogonal form
Nsites = length(A);
L{1} = ML; R{Nsites} = MR;
chid = size(M{1},3); 
for p = 1:Nsites - 1
    chil = size(A{p},1); chir = size(A{p},3);
    [qtemp,rtemp] = qr(reshape(A{p},[chil*chid,chir]),0);
    A{p} = reshape(qtemp,[chil,chid,chir]);
    A{p+1} = ncon({rtemp,A{p+1}},{[-1,1],[1,-2,-3]})/norm(rtemp(:));
    L{p+1} = ncon({L{p},M{p},A{p},conj(A{p})},{[2,1,4],[2,-1,3,5],[4,5,-3],[1,3,-2]});
end
chil = size(A{Nsites},1); chir = size(A{Nsites},3);
[qtemp,stemp] = qr(reshape(A{Nsites},[chil*chid,chir]),0);
A{Nsites} = reshape(qtemp,[chil,chid,chir]);
sWeight{Nsites+1} = stemp./sqrt(trace(stemp*stemp'));

Ekeep = [];
error=0;chimax=Band.chimin;
for k = 1:OPTS.numsweeps
    if error>1e-9
        chimax=min(chimax+Band.chi_step,Band.chimax);
    elseif error<1e-9
        chimax=max(chimax-Band.chi_step,Band.chimin);
    end
    error=0;
    Mode_A=1;
    %%%%%% Optimization sweep: right-to-left 
    for p = Nsites-1:-1:1
        
        %%%%% two-site update
        chil = size(A{p},1); chir = size(A{p+1},3);
        psiGround = reshape(ncon({A{p},A{p+1},sWeight{p+2}},{[-1,-2,1],[1,-3,2],[2,-4]}),[chil*chid^2*chir,1]);
        if OPTS.updateon==1
            [psiGround,Ekeep(end+1)] = eig_Arnoldi(psiGround,OPTS,@doApplyMPO,{L{p},M{p},M{p+1},R{p+1}});
        elseif OPTS.updateon==2
            [psiGround,Ekeep(end+1)] = eigLanczos(psiGround,OPTS,@doApplyMPO,{L{p},M{p},M{p+1},R{p+1}});
        end
        [utemp,stemp,vtemp] = svd(reshape(psiGround,[chil*chid,chid*chir]),'econ');
        chitemp = min(min(size(stemp)),chimax);
        A{p} = reshape(utemp(:,1:chitemp),[chil,chid,chitemp]);
        Mode=sqrt(sum(diag(stemp(1:chitemp,1:chitemp)).^2));Mode_all=sqrt(sum(diag(stemp).^2));
        error=error+(Mode_all-Mode)/Mode_all;
        Mode_A=Mode_A*Mode;
        sWeight{p+1} = stemp(1:chitemp,1:chitemp)./Mode;
        B{p+1} = reshape(vtemp(:,1:chitemp)',[chitemp,chid,chir]);
            
        %%%%% new block Hamiltonian MPO
        R{p} = ncon({M{p+1},R{p+1},B{p+1},conj(B{p+1})},{[-1,2,3,5],[2,1,4],[-3,5,4],[-2,3,1]});
        
       %%%%% display energy
        if OPTS.display == 2
            fprintf('Sweep: %2.1d of %2.1d, Loc: %2.1d, Energy: %12.12d\n',k,OPTS.numsweeps,p,Ekeep(end));
        end 
    end
    
    %%%%%% left boundary tensor
    chil = size(A{1},1); chir = size(A{1},3);
    [utemp,stemp,vtemp] = svd(reshape(ncon({A{1},sWeight{2}},{[-1,-2,1],[1,-3]}),[chil,chid*chir]),'econ');
    B{1} = reshape(vtemp',[chil,chid,chir]);
    Mode=sqrt(trace(stemp.^2));
    Mode_A=Mode_A*Mode;
    sWeight{1} = utemp*stemp./Mode;
    
    %%%%%% Optimization sweep: left-to-right
    for p = 1:Nsites-1
        
        %%%%% two-site update
        chil = size(B{p},1); chir = size(B{p+1},3);
        psiGround = reshape(ncon({sWeight{p},B{p},B{p+1}},{[-1,1],[1,-2,2],[2,-3,-4]}),[chil*chid^2*chir,1]);
        if OPTS.updateon==1 
            [psiGround,Ekeep(end+1)] = eig_Arnoldi(psiGround,OPTS,@doApplyMPO,{L{p},M{p},M{p+1},R{p+1}});
        elseif OPTS.updateon==2 
            [psiGround,Ekeep(end+1)] = eigLanczos(psiGround,OPTS,@doApplyMPO,{L{p},M{p},M{p+1},R{p+1}});
        end
        [utemp,stemp,vtemp] = svd(reshape(psiGround,[chil*chid,chid*chir]),'econ');
        chitemp = min(min(size(stemp)),chimax);
        A{p} = reshape(utemp(:,1:chitemp),[chil,chid,chitemp]);
        Mode=sqrt(sum(diag(stemp(1:chitemp,1:chitemp)).^2));Mode_all=sqrt(sum(diag(stemp).^2));
        error=error+(Mode_all-Mode)/Mode_all;
        Mode_A=Mode_A*Mode;
        sWeight{p+1} = stemp(1:chitemp,1:chitemp)./Mode;
        B{p+1} = reshape(vtemp(:,1:chitemp)',[chitemp,chid,chir]);
            
        %%%%% new block Hamiltonian
        L{p+1} = ncon({L{p},M{p},A{p},conj(A{p})},{[2,1,4],[2,-1,3,5],[4,5,-3],[1,3,-2]});
        
        %%%%% display energy
        if OPTS.display == 2
            fprintf('Sweep: %2.1d of %2.1d, Loc: %2.1d, Energy: %12.12d\n',k,OPTS.numsweeps,p,Ekeep(end));
        end
    end
    
    %%%%%% right boundary tensor
    chil = size(B{Nsites},1); chir = size(B{Nsites},3);
    [utemp,stemp,vtemp] = svd(reshape(ncon({B{Nsites},sWeight{Nsites}},{[1,-2,-3],[-1,1]}),[chil*chid,chir,1]),'econ');    
    A{Nsites} = reshape(utemp,[chil,chid,chir]);
    Mode=sqrt(trace(stemp.^2));
    Mode_A=Mode_A*Mode;
    sWeight{Nsites+1} = (stemp./Mode)*vtemp';
    
    if OPTS.display == 1
        fprintf('Sweep: %2.1d of %2.1d, CutEr:%2.2d, Mode:%2.2d, Energy: %12.12d, Bond dim: %2.0d\n',k,OPTS.numsweeps,error,Mode_A,Ekeep(end),chimax);
    end
end
% A{Nsites} = ncon({A{Nsites},sWeight{Nsites+1}},{[-1,-2,1],[1,-3]});
% sWeight{Nsites+1} = eye(size(A{Nsites},3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = doApplyMPO(psi,L,M1,M2,R)
% applies the superblock MPO to the state

psi = reshape(ncon({reshape(psi,[size(L,3),size(M1,4),size(M2,4),size(R,3)]),L,M1,M2,R},...
    {[1,3,5,7],[2,-1,1],[2,4,-2,3],[4,6,-3,5],[6,-4,7]}),[size(L,3)*size(M1,4)*size(M2,4)*size(R,3),1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psivec,dval] = eigLanczos(psivec,OPTS,linFunct,functArgs)
% function for computing the smallest algebraic eigenvalue and eigenvector
% of the linear function 'linFunct' using a Lanczos method. Maximum
% iterations are specified by 'OPTS.maxit' and the dimension of Krylov
% space is specified by 'OPTS.krydim'. Input 'functArgs' is an array of
% optional arguments passed to 'linFunct'.

if norm(psivec) == 0
    psivec = rand(length(psivec),1);
end
psi = zeros(numel(psivec),OPTS.krydim+1);
A = zeros(OPTS.krydim,OPTS.krydim);
for k = 1:OPTS.maxit
    
    psi(:,1) = psivec(:)/norm(psivec);
    for p = 2:OPTS.krydim+1
        psi(:,p) = linFunct(psi(:,p-1),functArgs{(1:length(functArgs))});
        for g = 1:1:p-1
            A(p-1,g) = dot(psi(:,p),psi(:,g));
            A(g,p-1) = conj(A(p-1,g));
        end
        for g = 1:1:p-1
            psi(:,p) = psi(:,p) - dot(psi(:,g),psi(:,p))*psi(:,g);
            psi(:,p) = psi(:,p)/max(norm(psi(:,p)),1e-16);
        end
    end
    
    [utemp,dtemp] = eig(A);
    Evalue=real(diag(dtemp));
    xloc = find(Evalue == min(Evalue));
    psivec = psi(:,1:OPTS.krydim)*utemp(:,xloc(1));
end
psivec = psivec/norm(psivec);
dval = dtemp(xloc(1),xloc(1));

function [psivec,dval] = eig_Arnoldi(psivec,OPTS,linFunct,functArgs)
% function for computing the smallest algebraic eigenvalue and eigenvector
% of the linear function 'linFunct' using a Lanczos method. Maximum
% iterations are specified by 'OPTS.maxit' and the dimension of Krylov
% space is specified by 'OPTS.krydim'. Input 'functArgs' is an array of
% 114514
% optional arguments passed to 'linFunct'.

if norm(psivec) == 0
    psivec = rand(length(psivec),1);
end
psi = zeros(numel(psivec),OPTS.krydim+1);
A = zeros(OPTS.krydim+1,OPTS.krydim);
for k = 1:OPTS.maxit
    psi(:,1) = psivec(:)/norm(psivec);
    for p = 2:OPTS.krydim+1
        psi0 = linFunct(psi(:,p-1),functArgs{(1:length(functArgs))});
        psi(:,p) = psi0;
        for g = 1:1:p-1
            A(g,p-1)=dot(psi(:,g),psi0);
        end
        for g = 1:1:p-1
            h=dot(psi(:,g),psi(:,p));
            psi(:,p) = psi(:,p) - h*psi(:,g);
            psi(:,p) = psi(:,p)/max(norm(psi(:,p)),1e-16);
        end
        A(p,p-1)=dot(psi(:,p),psi0);
    end

    [Rutemp,dtemp,~] = eig(A(1:OPTS.krydim,1:OPTS.krydim));
    Evalue=real(diag(dtemp));
    xloc = find(Evalue == min(Evalue));
    psivec = psi(:,1:OPTS.krydim)*Rutemp(:,xloc(1));
end
psivec = psivec/norm(psivec);
dval = dtemp(xloc(1),xloc(1));