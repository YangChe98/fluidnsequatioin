function [erroru,errorp]=H1_semierror(x_n,y_n,t_ini,tn,deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana,basisfunctionweighttotal,node_coordinate0,element_coordinate0,node_coordinate,element_coordinate,localbasisfunctionnumber1,localbasisfunctionnumber2,ubasis_function_number,node_number)
[element_number,n]=size(element_coordinate0);
load phigaussvaluedx.mat;
load phigaussvaluedy.mat;
load phigaussvaluedxorder1.mat;
load phigaussvaluedyorder1.mat;

gaussweight=[5/9,8/9,5/9];
gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

xgausspoint2d=repmat(gausspoint,3,1);
ygausspoint2d=repmat(gausspoint.',1,3);
xgausspoint2d=reshape(xgausspoint2d,1,[]);
ygausspoint2d=reshape(ygausspoint2d,1,[]);

gaussweight2d=gaussweight.*gaussweight.';
gaussweight2d=reshape(gaussweight2d,1,[]);
gaussweight2d1=repmat(gaussweight2d,localbasisfunctionnumber1,1);
gaussweight2d2=repmat(gaussweight2d,localbasisfunctionnumber2,1);

J=[1/(x_n*2),0;0,1/(y_n*2)];
 detJ=det(J);
 detJ=abs(detJ);
 invJ=inv(J);
 invJT=invJ.';
%  C=invJ*invJT;
 phigaussvaluedxlocal= phigaussvaluedx*J(1,1)+ phigaussvaluedy*J(2,1);
  phigaussvaluedylocal= phigaussvaluedx*J(1,2)+ phigaussvaluedy*J(2,2);
   phigaussvaluedxlocalorder1= phigaussvaluedxorder1*J(1,1)+ phigaussvaluedyorder1*J(2,1);
  phigaussvaluedylocalorder1= phigaussvaluedxorder1*J(1,2)+ phigaussvaluedyorder1*J(2,2);
erroru=0;
errorp=0;
t_now=t_ini;

for t=1:tn

    
for i=1:element_number
     local_coordinate1=node_coordinate(element_coordinate(i,1),1:2);
%     local_coordinate2=node_coordinate(element_coordinate(i,2),1:2);
    local_coordinate3=node_coordinate(element_coordinate(i,3),1:2);
%     local_coordinate4=node_coordinate(element_coordinate(i,4),1:2);
    
    Fkb=(local_coordinate1+local_coordinate3).'/2;
   
    xlocalgausspoint=J(1,1)*xgausspoint2d+J(1,2)*ygausspoint2d+Fkb(1,1)*ones(size(xgausspoint2d));
    ylocalgausspoint=J(2,1)*xgausspoint2d+J(2,2)*ygausspoint2d+Fkb(2,1)*ones(size(xgausspoint2d));
  weightlocal1=full(basisfunctionweighttotal(element_coordinate(i,:),t));
    weightlocal2=full(basisfunctionweighttotal(element_coordinate(i,:)+node_number,t));
    weightlocal3=full(basisfunctionweighttotal(element_coordinate0(i,:)+ubasis_function_number,t));
%     ux=f(xlocalgausspoint,ylocalgausspoint,t_now,nu);
    [ux,uy,vx,vy,px,py]=function_uvpderivative(xlocalgausspoint,ylocalgausspoint,t_now);
    
    
    phixlocal1=weightlocal1.'*phigaussvaluedxlocal;
    phiylocal1=weightlocal1.'*phigaussvaluedylocal;
    phixlocal2=weightlocal2.'*phigaussvaluedxlocal;
    phiylocal2=weightlocal2.'*phigaussvaluedylocal;
     phixlocal3=weightlocal3.'*phigaussvaluedxlocalorder1;
    phiylocal3=weightlocal3.'*phigaussvaluedylocalorder1;
   
    
    errorp=errorp+deltat*detJ*sqrt((px-phixlocal3).^2*gaussweight2d.')+deltat*detJ*sqrt((py-phiylocal3).^2*gaussweight2d.');
    erroru=erroru+deltat*detJ*sqrt((ux-phixlocal1).^2*gaussweight2d.')+deltat*detJ*sqrt((uy-phiylocal1).^2*gaussweight2d.')...
        +deltat*detJ*sqrt((vx-phixlocal2).^2*gaussweight2d.')+deltat*detJ*sqrt((vy-phiylocal2).^2*gaussweight2d.');
end
end
t_now=t_now+deltat;
end