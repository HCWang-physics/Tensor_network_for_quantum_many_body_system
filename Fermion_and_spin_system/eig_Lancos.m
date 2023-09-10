function [psivec,dval,A] = eig_Lancos(psivec,H,OPTS)
% function for computing the smallest algebraic eigenvalue and eigenvector
% of the linear function 'linFunct' using a Lanczos method. Maximum
% iterations are specified by 'OPTS.maxit' and the dimension of Krylov
% space is specified by 'OPTS.krydim'. Input 'functArgs' is an array of
% optional arguments passed to 'linFunct'.
dval=zeros(1,OPTS.maxit);

if norm(psivec) == 0
    psivec = rand(length(psivec),1);
end
psi = zeros(numel(psivec),OPTS.krydim+1);
A = zeros(OPTS.krydim,OPTS.krydim);
for k = 1:OPTS.maxit
    psi(:,1) = psivec(:)/norm(psivec);
    for p = 2:OPTS.krydim+1
        psi(:,p) = H*psi(:,p-1);
        for g = 1:1:p-1
            A(p-1,g) = dot(psi(:,p),psi(:,g));
            A(g,p-1) = conj(A(p-1,g));
        end
        for g = 1:1:p-1
            psi(:,p) = psi(:,p) - dot(psi(:,g),psi(:,p))*psi(:,g);
            psi(:,p) = psi(:,p)/max(norm(psi(:,p)),1e-16);
        end
    end

    [Rutemp,dtemp] = eig(A);
    Evalue=real(diag(dtemp));
    xloc = find(Evalue == min(Evalue));
    psivec = psi(:,1:OPTS.krydim)*Rutemp(:,xloc(1));
    dval(1,k) = dtemp(xloc(1),xloc(1));
end
psivec = psivec/norm(psivec);
end