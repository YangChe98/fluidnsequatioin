% function []=fexpression()

syms x y t nu

u1=(x.^2.*y.^2+exp(-y)).*cos(2*pi*t);
u2=(-2*x.*y.^3/3+2-pi*sin(pi*x)).*cos(2*pi*t);
p=-(2-pi*sin(pi*x)).*cos(2*pi*y).*cos(2*pi*t);
% u1=t + x.^2 + y.^2;
% u2=t + 2*x.^2 - 2*x.*y;
% p= x + y - 1;
u1t=diff(u1,t);
u2t=diff(u2,t);
u1x=diff(u1,x);
u2x=diff(u2,x);
u1y=diff(u1,y);
u2y=diff(u2,y);
u1xx=diff(u1,x,2);
u2xx=diff(u2,x,2);
u1yy=diff(u1,y,2);
u2yy=diff(u2,y,2);
px=diff(p,x);
py=diff(p,y);
f1=u1t+(u1*u1x+u2*u1y)+px-nu*(u1xx+u1yy);
f2=u2t+(u1*u2x+u2*u2y)+py-nu*(u2xx+u2yy);