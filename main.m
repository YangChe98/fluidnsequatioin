
%% 2D Quadrilateral Finite Element Methods Using Discontinuous Pressures
% u velocity using Q2 pressure using Q1 ,[X,Y] is the mesh of u, [X0, Y0]
% is the mesh of p. 0 is for p. without 0 is for u.

clc
clear
tic
nu=1;
syms xi eta 


x_max=1;
x_min=0;
y_max=1;
y_min=0;
x_n=2^4;
y_n=2^4;
t_ini=0;
t_end=1;
deltat=1e-2;
deltax=(x_max-x_min)/x_n;
deltay=(y_max-y_min)/y_n;
t=t_ini:deltat:t_end;
tn=floor((t_end-t_ini)/deltat);

localbasisfunctionnumber1=9;
localbasisfunctionnumber2=4;
node_number=(2*x_n+1)*(2*y_n+1);
element_number=x_n*y_n;
ubasis_function_number=node_number*2;
pbasis_function_number=(x_n+1)*(y_n+1);
t0=0;



[X0,Y0,node_coordinate0,element_coordinate0,X,Y,node_coordinate,element_coordinate]=meshgenerate(x_min,x_max,y_min,y_max,x_n,y_n,t_ini);
pbasisfunctionweight0=function_p(node_coordinate0(:,1),node_coordinate0(:,2),t0);
u1basisfunctionweight0=function_u(node_coordinate(:,1),node_coordinate(:,2),t0);
u2basisfunctionweight0=function_v(node_coordinate(:,1),node_coordinate(:,2),t0);
basisfunctionweight0=[u1basisfunctionweight0;u2basisfunctionweight0;pbasisfunctionweight0];
basisfunctionweight0=sparse(basisfunctionweight0);


basisfunctionweighttotal=solution(deltat,t_ini,t_end,basisfunctionweight0,nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number,x_min,x_max,y_min,y_max);
uvalue=basisfunctionweighttotal(1:node_number,:);

vvalue=basisfunctionweighttotal(node_number+1:node_number*2,:);


pvalue=basisfunctionweighttotal(node_number*2+1:ubasis_function_number+pbasis_function_number,:);

uvalueana=function_u(node_coordinate(:,1),node_coordinate(:,2),t);
vvalueana=function_v(node_coordinate(:,1),node_coordinate(:,2),t);
pvalueana=function_p(node_coordinate0(:,1),node_coordinate0(:,2),t);
plotuvps(tn,x_min,x_n,x_max,y_min,y_n,y_max,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana)
%  plotuvp0(x_min,x_n,x_max,y_min,y_n,y_max,basisfunctionweight,node_number,ubasis_function_number,pbasis_function_number)
[errorLinfity_u,errorLinfity_p]=Linfinity_error(deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana);
[errorL2_u,errorL2_p,errorH1_u,errorH1_p]=H1_error(x_n,y_n,t_ini,tn,deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana,basisfunctionweighttotal,node_coordinate0,element_coordinate0,node_coordinate,element_coordinate,localbasisfunctionnumber1,localbasisfunctionnumber2,ubasis_function_number,node_number);

toc