function [phigaussvaluedxorder1,phigaussvaluedyorder1]=basisgausspointvaluederivativeorder1()

syms x y;

gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

xgausspoint2d=repmat(gausspoint,3,1);
ygausspoint2d=repmat(gausspoint.',1,3);
xgausspoint2d=reshape(xgausspoint2d,1,[]);
ygausspoint2d=reshape(ygausspoint2d,1,[]);


  phigaussvaluedxorder1=zeros(4,9);

  phigaussvaluedyorder1=zeros(4,9);
for i=1:4
  
     [phi,phix,phiy]=Q1basisfunction(i,xgausspoint2d,ygausspoint2d);
 phigaussvaluedxorder1(i,:)=phix;
  phigaussvaluedyorder1(i,:)=phiy;

end

save phigaussvaluedxorder1.mat phigaussvaluedxorder1;
save phigaussvaluedyorder1.mat phigaussvaluedyorder1;