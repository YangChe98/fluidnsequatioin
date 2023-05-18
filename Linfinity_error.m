function [erroru,errorp]=Linfinity_error(deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana)
[m,n]=size(uvalue);
erroru=0;
errorp=0;
for i=1:n
    errorp=errorp+deltat*((norm(pvalueana(:,i)-pvalue(:,i),"inf"))*deltax*deltay);
    erroru=erroru+deltat*(max(abs(uvalueana(:,i)-uvalue(:,i))+abs(vvalueana(:,i)-vvalue(:,i))))*deltax*deltay/4;
end