function [ux,uy,vx,vy,px,py]=function_uvpderivative(x,y,t)
ux=2*x.*y.^2.*cos(2*pi*t);
vx=-cos(2*pi*t).*((2*y.^3)/3 + pi^2*cos(x*pi));
uy=-cos(2*pi*t).*(- 2*y.*x.^2 + exp(-y));
vy=-2*x.*y.^2.*cos(2*pi*t);
px=pi^2*cos(2*pi*t).*cos(pi*x).*cos(2*pi*y);
py=-2*pi*cos(2*pi*t).*sin(2*pi*y).*(pi*sin(pi*x) - 2);
