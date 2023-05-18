function [phiiphijphikxreference,phiiphijphikyreference]=referencematrixphiiphijphiky()

syms xi eta;

phiiphijphikxreference=zeros(9,9,9);
phiiphijphikyreference=zeros(9,9,9);
%A12reference=zeros(9,9);
%A22reference=zeros(9,9);


% reference matrix A 
for i=1:9
       phii=Q2basisfunction(i,xi,eta);
    for j=1:9
     
        phij=Q2basisfunction(j,xi,eta);
        for k=1:9
            [phik,phikx,phiky]=Q2basisfunction(k,xi,eta);
        A11=int(matlabFunction(phii.*phij.*phikx,'vars',[xi eta]),xi,-1,1);
        A1=int(matlabFunction(A11,'var',[eta]),eta,-1,1);
        phiiphijphikxreference(i,j,k)=double(A1);
%         


          A22=int(matlabFunction(phii.*phij.*phiky,'vars',[xi eta]),xi,-1,1);
        A2=int(matlabFunction(A22,'var',[eta]),eta,-1,1);
        phiiphijphikyreference(i,j,k)=double(A2);

        end
    end
end

save phiiphijphikxreference.mat phiiphijphikxreference;
save phiiphijphikyreference.mat phiiphijphikyreference;
end