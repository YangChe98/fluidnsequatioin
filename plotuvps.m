function plotuvps(tn,x_min,x_n,x_max,y_min,y_n,y_max,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana)

xplot=x_min:(1/(2*x_n)):x_max;
 yplot=y_min:(1/(2*y_n)):y_max;
xplot1=x_min:(1/(x_n)):(x_max);
 yplot1=y_min:(1/(y_n)):(y_max);

[xplot,yplot]=meshgrid(xplot,yplot);
[xplot1,yplot1]=meshgrid(xplot1,yplot1);

for i=1:(tn/50+1)
uvalueplot=uvalue(:,1+(i-1)*50);
uvalueplot=reshape(uvalueplot,size(xplot.'));
vvalueplot=vvalue(:,1+(i-1)*50);
vvalueplot=reshape(vvalueplot,size(xplot.'));

pvalueplot=pvalue(:,1+(i-1)*50);
 pvalueplot=reshape(pvalueplot,size(xplot1.'));

 uvalueanaplot=uvalueana(:,1+(i-1)*50);
uvalueanaplot=reshape(uvalueanaplot,size(xplot.'));
vvalueanaplot=vvalueana(:,1+(i-1)*50);
vvalueanaplot=reshape(vvalueanaplot,size(xplot.'));

pvalueanaplot=pvalueana(:,1+(i-1)*50);
 pvalueanaplot=reshape(pvalueanaplot,size(xplot1.'));
figure(1+(i-1)*8)
contourf(xplot,yplot,sqrt(uvalueplot.'.^2+vvalueplot.'.^2))
l1=streamslice(xplot,yplot,uvalueplot.',vvalueplot.');
set(l1,'Color','w');

figure(2+(i-1)*8)

mesh(xplot1,yplot1,pvalueplot.')
colorbar
figure(3+(i-1)*8)

mesh(xplot,yplot,uvalueplot.')
colorbar
figure(4+(i-1)*8)

mesh(xplot,yplot,vvalueplot.')
colorbar

figure(5+(i-1)*8)

mesh(xplot1,yplot1,pvalueanaplot.')
colorbar
figure(6+(i-1)*8)

mesh(xplot,yplot,uvalueanaplot.')
colorbar
figure(7+(i-1)*8)

mesh(xplot,yplot,vvalueanaplot.')
figure(8+(i-1)*8)
contourf(xplot,yplot,sqrt(uvalueanaplot.'.^2+vvalueanaplot.'.^2))
l2=streamslice(xplot,yplot,uvalueanaplot.',vvalueanaplot.');
set(l2,'Color','w');
end