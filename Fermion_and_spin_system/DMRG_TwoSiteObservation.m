function [Twosite_correlation] = DMRG_TwoSiteObservation(F,A,sWeight,OperS,OperI,index,Sele_site)
%UNTITLED8 此处提供此函数的摘要
%   此处提供详细说明
Site_number=length(A);
Twosite_correlation=zeros(1,Site_number);
if index==1 %spin up
    OperS1=OperS;OperI1=OperI*F;
    OperS2=-OperS*F;OperI2=OperI;
    %main progress
    Psi=ncon({A{Sele_site},OperS1,conj(A{Sele_site}),sWeight{Sele_site+1},conj(sWeight{Sele_site+1})},...
        {[-1,4,2],[5,4],[-2,5,3],[2,1],[3,1]});
    for j=Sele_site-1:-1:1
        Twosite_correlation(j)=ncon({Psi,A{j},OperI1,conj(A{j})},{[1,4],[5,2,1],[3,2],[5,3,4]});
        Psi=ncon({Psi,A{j},F,conj(A{j})},{[1,4],[-1,2,1],[3,2],[-2,3,4]});
    end
    
    Twosite_correlation(Sele_site)=ncon({A{Sele_site},OperI*OperS,conj(A{Sele_site}), ...
        sWeight{Sele_site+1},conj(sWeight{Sele_site+1})},{[1,2,4],[3,2],[1,3,5],[4,6],[5,6]});
    
    Psi=ncon({A{Sele_site},OperS2,conj(A{Sele_site})},{[1,2,-1],[3,2],[1,3,-2]});
    for j=Sele_site+1:Site_number
        Twosite_correlation(j)=ncon({Psi,A{j},OperI2,conj(A{j}),sWeight{j+1},conj(sWeight{j+1})},{[1,4],[1,2,5],[3,2],[4,3,6],[5,7],[6,7]});%{[4,5],[4,6,2],[7,6],[5,7,3],[2,1],[3,1]}
        Psi=ncon({Psi,A{j},F,conj(A{j})},{[1,4],[1,2,-1],[3,2],[4,3,-2]});
    end

elseif index==2 %spin correlation
    Psi=ncon({A{Sele_site},OperS,conj(A{Sele_site}),sWeight{Sele_site+1},conj(sWeight{Sele_site+1})},...
        {[-1,4,2],[5,4],[-2,5,3],[2,1],[3,1]});
    for j=Sele_site-1:-1:1
        Twosite_correlation(j)=ncon({Psi,A{j},OperI,conj(A{j})},{[1,4],[5,2,1],[3,2],[5,3,4]});
        Psi=ncon({Psi,A{j},conj(A{j})},{[1,3],[-1,2,1],[-2,2,3]});
    end
    
    Twosite_correlation(Sele_site)=ncon({A{Sele_site},OperI*OperS,conj(A{Sele_site}), ...
        sWeight{Sele_site+1},conj(sWeight{Sele_site+1})},{[1,2,4],[3,2],[1,3,5],[4,6],[5,6]});
    
    Psi=ncon({A{Sele_site},OperS,conj(A{Sele_site})},{[1,2,-1],[3,2],[1,3,-2]});
    for j=Sele_site+1:Site_number
        Twosite_correlation(j)=ncon({Psi,A{j},OperI,conj(A{j}),sWeight{j+1},conj(sWeight{j+1})},{[1,4],[1,2,5],[3,2],[4,3,6],[5,7],[6,7]});
        Psi=ncon({Psi,A{j},conj(A{j})},{[1,3],[1,2,-1],[3,2,-2]});
    end
end

end