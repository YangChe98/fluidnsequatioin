function [g1,g2]=boundarygt(t,x,y)
% g1=t + x.^2 + y.^2;
% g2=t + 2*x.^2 - 2*x.*y;
g1=function_u(x,y,t);
g2=function_v(x,y,t);
end