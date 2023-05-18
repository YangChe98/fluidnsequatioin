function vvalue=function_v(x,y,t)
%vvalue=-pi*cos(pi*x).*y;
%vvalue= -2/3*x.*y.^3+2-pi*sin(pi*x);
vvalue=(-2*x.*y.^3/3+2-pi*sin(pi*x)).*cos(2*pi*t);
 %vvalue=t + 2*x.^2 - 2*x.*y;
end