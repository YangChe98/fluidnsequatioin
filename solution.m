function basisfunctionweighttotal=solution(deltat,t_ini,t_end,basisfunctionweight0,nu,x_n,y_n,element_number,localbasisfunctionnumber1,element_coordinate0,element_coordinate,node_coordinate,ubasis_function_number,pbasis_function_number,node_number,x_min,x_max,y_min,y_max)


timesteptotal=floor((t_end-t_ini)/deltat);
basisfunctionweighttotal=zeros(ubasis_function_number+pbasis_function_number,timesteptotal+1);
basisfunctionweightprevious=basisfunctionweight0;
basisfunctionweighttotal(:,1)=basisfunctionweight0;

load A11reference.mat;
load A12reference.mat;
load A22reference.mat;

load B11reference.mat;
load B12reference.mat;



load phiiphijphikxreference.mat;
load phiiphijphikyreference.mat
load phiiphijreference.mat;

load phigaussvalue.mat;


%%%%%%%%%%%%%%% gauss quadrature %%%%%%%%%%
gaussweight=[5/9,8/9,5/9];
gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

xgausspoint2d=repmat(gausspoint,3,1);
ygausspoint2d=repmat(gausspoint.',1,3);
xgausspoint2d=reshape(xgausspoint2d,1,[]);
ygausspoint2d=reshape(ygausspoint2d,1,[]);

gaussweight2d=gaussweight.*gaussweight.';
gaussweight2d=reshape(gaussweight2d,1,[]);
gaussweight2d=repmat(gaussweight2d,localbasisfunctionnumber1,1);

A=sparse(ubasis_function_number,ubasis_function_number);
B=sparse(ubasis_function_number,pbasis_function_number);
right=sparse(ubasis_function_number+pbasis_function_number,1);
imatrix=[];


J=[1/(x_n*2),0;0,1/(y_n*2)];
 detJ=det(J);
 detJ=abs(detJ);
 invJ=inv(J);
 invJT=invJ.';
 C=invJ*invJT;
 Alocal=nu*detJ*(C(1,1)*(A11reference)+C(1,2)*(A12reference+A12reference.')+C(2,2)*(A22reference));
    B1local=detJ*(invJ(1,1)*B11reference+invJ(2,1)*B12reference);
      B2local=detJ*(invJ(1,2)*B11reference+invJ(2,2)*B12reference);
      phiiphij=detJ*phiiphijreference;
     phiiphijphikx=detJ*( phiiphijphikxreference*invJ(1,1)+phiiphijphikyreference*invJ(2,1));
      phiiphijphiky= detJ*(phiiphijphikxreference*invJ(1,2)+phiiphijphikyreference*invJ(2,2));
%%%%%%%%%%%%gengerate  matrix A   B right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_now=t_ini;

for tn=1:timesteptotal

    A=sparse(ubasis_function_number,ubasis_function_number);
B=sparse(ubasis_function_number,pbasis_function_number);
right=sparse(ubasis_function_number+pbasis_function_number,1);


    t_now=t_now+deltat;
    node_coordinate=meshboundarygenerate(node_coordinate,t_now,x_min,x_max,y_min,y_max);
for i=1:element_number 
    phiiphijphikxsumkweight1=zeros(9,9);
       phiiphijphikxsumkweight2=zeros(9,9);
        phiiphijphikxsumiweight1=zeros(9,9);
       phiiphijphikxsumiweight2=zeros(9,9);
  phiiphijphikysumkweight1=zeros(9,9);
       phiiphijphikysumkweight2=zeros(9,9);
        phiiphijphikysumiweight1=zeros(9,9);
       phiiphijphikysumiweight2=zeros(9,9);

    local_coordinate1=node_coordinate(element_coordinate(i,1),1:2);
    local_coordinate2=node_coordinate(element_coordinate(i,2),1:2);
    local_coordinate3=node_coordinate(element_coordinate(i,3),1:2);
    local_coordinate4=node_coordinate(element_coordinate(i,4),1:2);
    
    Fkb=(local_coordinate1+local_coordinate3).'/2;
   
    xlocalgausspoint=J(1,1)*xgausspoint2d+J(1,2)*ygausspoint2d+Fkb(1,1)*ones(size(xgausspoint2d));
    ylocalgausspoint=J(2,1)*xgausspoint2d+J(2,2)*ygausspoint2d+Fkb(2,1)*ones(size(xgausspoint2d));
 
    fvalue=f(xlocalgausspoint,ylocalgausspoint,t_now,nu);

    
     
   
    rightlocal1=detJ*(fvalue(1,:).*phigaussvalue).*gaussweight2d;
    rightlocal2=detJ*(fvalue(2,:).*phigaussvalue).*gaussweight2d;
    rightlocal1=sum(rightlocal1,2);
    rightlocal2=sum(rightlocal2,2);
    weightlocal1=full(basisfunctionweightprevious(element_coordinate(i,:),:));
    weightlocal2=full(basisfunctionweightprevious(element_coordinate(i,:)+node_number,:));

 phiiphijlocal1=phiiphij*weightlocal1;
 phiiphijlocal2=phiiphij*weightlocal2;


 for k=1:9

      phiiphijphikxsumkweight1= phiiphijphikxsumkweight1+squeeze(phiiphijphikx(:,:,k)*weightlocal1(k));
     
%        phiiphijphikxsumkweight2= phiiphijphikxsumkweight2+squeeze(phiiphijphikx(:,:,k)*weightlocal2(k));
     phiiphijphikxsumiweight1= phiiphijphikxsumiweight1+squeeze(phiiphijphikx(k,:,:)*weightlocal1(k));
%   phiiphijphikxsumiweight2= phiiphijphikxsumiweight2+squeeze(phiiphijphikx(k,:,:)*weightlocal2(k));

   phiiphijphikysumkweight1= phiiphijphikysumkweight1+squeeze(phiiphijphiky(:,:,k)*weightlocal1(k));
     
       phiiphijphikysumkweight2= phiiphijphikysumkweight2+squeeze(phiiphijphiky(:,:,k)*weightlocal2(k));
     phiiphijphikysumiweight1= phiiphijphikysumiweight1+squeeze(phiiphijphiky(k,:,:)*weightlocal1(k));
  phiiphijphikysumiweight2= phiiphijphikysumiweight2+squeeze(phiiphijphiky(k,:,:)*weightlocal2(k));
 end

A(element_coordinate(i,:),element_coordinate(i,:))=A(element_coordinate(i,:),element_coordinate(i,:))+deltat*Alocal+phiiphij+deltat*phiiphijphikxsumiweight1+deltat*phiiphijphikysumiweight2+deltat*phiiphijphikxsumkweight1/2+deltat*phiiphijphikysumkweight1/2;
 A(element_coordinate(i,:)+node_number,element_coordinate(i,:)+node_number)=A(element_coordinate(i,:)+node_number,element_coordinate(i,:)+node_number)+deltat*Alocal+phiiphij+deltat*phiiphijphikxsumiweight1+deltat*phiiphijphikysumiweight2+deltat*phiiphijphikxsumkweight2/2+deltat*phiiphijphikysumkweight2/2;
 B(element_coordinate(i,:),element_coordinate0(i,:))=B(element_coordinate(i,:),element_coordinate0(i,:))+B1local;
 B(element_coordinate(i,:)+node_number,element_coordinate0(i,:))=B(element_coordinate(i,:)+node_number,element_coordinate0(i,:))+B2local;
  right(element_coordinate(i,:),1) =right(element_coordinate(i,:),1)+deltat*rightlocal1+phiiphijlocal1;
   right(element_coordinate(i,:)+node_number,1) =right(element_coordinate(i,:)+node_number,1)+deltat*rightlocal2+phiiphijlocal2;
end

B=deltat*B;
Atotal=[A,B;B.',zeros(pbasis_function_number,pbasis_function_number)];

%%%%%%%%%%%%%%%%%%%%% boundary condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:node_number
    if(node_coordinate(i,3)==1)
        right=right-Atotal(:,i)*node_coordinate(i,4)-Atotal(:,node_number+i)*node_coordinate(i,5);
        Atotal(:,i)=0;
        Atotal(i,:)=0;
        
        Atotal(:,node_number+i)=0;
        Atotal(node_number+i,:)=0;
        Atotal(i,i)=1;
        Atotal(node_number+i,node_number+i)=1;
        right(i,1)=node_coordinate(i,4);
        right(node_number+i,1)=node_coordinate(i,5);
    end

end

 right=right-Atotal(:,ubasis_function_number+1)*node_coordinate(1,6);
        Atotal(:,ubasis_function_number+1)=0;
        Atotal(ubasis_function_number+1,:)=0;
        Atotal(ubasis_function_number+1,ubasis_function_number+1)=1;
        right(ubasis_function_number+1,1)=node_coordinate(1,6);
 

%%%%%%%%%%%%%%%%%% solution %%%%%%%%%%%%%%
basisfunctionweight=Atotal\right;
basisfunctionweightprevious=basisfunctionweight;
basisfunctionweighttotal(:,tn+1)=basisfunctionweight;
end

