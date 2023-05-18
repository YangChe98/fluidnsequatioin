function [phiiphijreference]=referencematrixphiiphij()

syms xi eta;

phiiphijreference=zeros(9,9);
%A12reference=zeros(9,9);
%A22reference=zeros(9,9);


% reference matrix A 
for i=1:9
       phii=Q2basisfunction(i,xi,eta);
    for j=1:9
     
        phij=Q2basisfunction(j,xi,eta);
        A11=int(matlabFunction(phii.*phij,'vars',[xi eta]),xi,-1,1);
        A1=int(matlabFunction(A11,'var',[eta]),eta,-1,1);
        phiiphijreference(i,j)=double(A1);
%         


    end
end

save phiiphijreference.mat phiiphijreference;

end