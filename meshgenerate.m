function [X0,Y0,node_coordinate0,element_coordinate0,X,Y,node_coordinate,element_coordinate]=meshgenerate(x_min,x_max,y_min,y_max,x_n,y_n,t_ini)
% local coordinate Q2
% 4--7--3
% |  |  |
% 8--9--6
% |  |  |
% 1--5--2
% local coordinate Q1
% 4-- --3
% |  |  |
%  -- --
% |  |  |
% 1-- --2

%global coordinate: Q2
%21--22-23-24-25
%|      |     | 
%16  17 18 19 20
%|      |     |
%11--12-13-14-15
%|      |     |
%6   7  8  9  10
%|      |     |
% 1--2--3--4--5
%global coordinate: Q1
%7------8-----9
%|      |     | 
%|      |     |
%4------5-----6
%|      |     |
%|      |     |
% 1-----2-----3
deltax=(x_max-x_min)/x_n;
deltay=(y_max-y_min)/y_n;
x0=x_min:deltax:x_max;
y0=y_min:deltay:y_max;

[X0,Y0]=meshgrid(x0,y0);
X0=reshape(X0',[],1);
Y0=reshape(Y0',[],1);
node_coordinate0=[X0,Y0];


x=x_min:deltax/2:x_max;
y=y_min:deltay/2:y_max;
[X,Y]=meshgrid(x,y);
X=reshape(X',[],1);
Y=reshape(Y',[],1);
node_coordinate=[X,Y];
[m,n]=size(X);

node_coordinate=meshboundarygenerate(node_coordinate,t_ini,x_min,x_max,y_min,y_max);



element_coordinate11=(1:x_n)';

element_coordinate1=[];
for i=1:y_n
    element_coordinate1=[element_coordinate1;element_coordinate11];
    element_coordinate11=element_coordinate11+x_n+1;
 
end
element_coordinate0=[element_coordinate1,element_coordinate1+1,element_coordinate1+x_n+2,element_coordinate1+x_n+1];



element_coordinate11=(1:2:2*x_n)';
element_coordinate1=[];
for i=1:y_n
    element_coordinate1=[element_coordinate1;element_coordinate11];
    element_coordinate11=element_coordinate11+2*(2*x_n+1);
end
element_coordinate=[element_coordinate1,element_coordinate1+2,element_coordinate1+2*(2*x_n+1)+2,element_coordinate1+2*(2*x_n+1),...
   element_coordinate1+1,element_coordinate1+(2*x_n+1)+2, element_coordinate1+2*(2*x_n+1)+1,element_coordinate1+(2*x_n+1),element_coordinate1+(2*x_n+1)+1];


