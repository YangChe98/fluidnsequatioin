function pvalue=function_p(x,y,t)
%pvalue=sin(2*pi*x)+sin(2*pi*y);
% pvalue=-cos(2*pi*y)*(2-pi*sin(pi*x));
pvalue=-(2-pi*sin(pi*x)).*cos(2*pi*y).*cos(2*pi*t);
%pvalue=( x + y - 1)*ones(size(t));
end