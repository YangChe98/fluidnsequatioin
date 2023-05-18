function uvalue=function_u(x,y,t)
%uvalue=sin(pi*x)+sin(pi*y);
% uvalue=x.^2.*y.^2 + exp(-y);
uvalue=(x.^2.*y.^2+exp(-y)).*cos(2*pi*t);
%uvalue=t + x.^2 + y.^2;
end