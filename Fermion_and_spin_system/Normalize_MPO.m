function [M] = Normalize_MPO(M)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
Nsites=length(M);
for p=1:Nsites-1
    Mtilda = permute(M{p},[4,1,3,2]);%change the index order of M
    chi1 = size(Mtilda,1); chi2 = size(Mtilda,2); chi3 = size(Mtilda,3); chi4 = size(Mtilda,4);
    [qtemp,rtemp] = qr(reshape(Mtilda,[chi1*chi2*chi3,chi4]),0);
    scaler=norm(rtemp(:));
    Mtilda=reshape(qtemp*scaler,[chi1,chi2,chi3,chi4]);
    M{p} = permute(Mtilda,[2,4,3,1]);%change the index order of M back to inital one
    M{p+1}=ncon({rtemp/scaler,M{p+1}},{[-1,1],[1,-2,-3,-4]});
end

for p=Nsites:2
    Mtilda = permute(M{p},[2,4,3,1]);%change the index order of M
    chi1 = size(Mtilda,1); chi2 = size(Mtilda,2); chi3 = size(Mtilda,3); chi4 = size(Mtilda,4);
    [qtemp,rtemp] = qr(reshape(Mtilda,[chi1*chi2*chi3,chi4]),0);
    scaler=norm(rtemp(:));
    Mtilda=reshape(qtemp*scaler,[chi1,chi2,chi3,chi4]);
    M{p} = permute(Mtilda,[4,1,3,2]);%change the index order of M back to inital one
    M{p-1}=ncon({rtemp/scaler,M{p-1}},{[1,-2],[-1,1,-3,-4]});
end
end