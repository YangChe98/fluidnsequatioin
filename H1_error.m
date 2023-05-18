function [errorL2_u,errorL2_p,errorH1_u,errorH1_p]=H1_error(x_n,y_n,t_ini,tn,deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana,basisfunctionweighttotal,node_coordinate0,element_coordinate0,node_coordinate,element_coordinate,localbasisfunctionnumber1,localbasisfunctionnumber2,ubasis_function_number,node_number)

[errorL2_u,errorL2_p]=L2_error(deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana);
[errorH1semi_u,errorH1semi_p]=H1_semierror(x_n,y_n,t_ini,tn,deltat,deltax,deltay,uvalue,vvalue,pvalue,uvalueana,vvalueana,pvalueana,basisfunctionweighttotal,node_coordinate0,element_coordinate0,node_coordinate,element_coordinate,localbasisfunctionnumber1,localbasisfunctionnumber2,ubasis_function_number,node_number);
errorH1_u=sqrt(errorL2_u^2+errorH1semi_u^2);
errorH1_p=sqrt(errorL2_p^2+errorH1semi_p^2);