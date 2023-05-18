function [erroru,errorp]=L2_error(deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana)
[m,n]=size(uvalue);
erroru=0;
errorp=0;
for i=1:n
    errorp=errorp+deltat*sqrt((norm(pvalueana(:,i)-pvalue(:,i))^2)*deltax*deltay);
    erroru=erroru+deltat*sqrt((norm(uvalueana(:,i)-uvalue(:,i))^2+norm(vvalueana(:,i)-vvalue(:,i))^2)*deltax*deltay/4);
end