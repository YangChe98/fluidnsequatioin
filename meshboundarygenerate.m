function node_coordinate=meshboundarygenerate(node_coordinate,t,x_min,x_max,y_min,y_max)
[m,n]=size(node_coordinate);
for i=1:m
    if(node_coordinate(i,2)==y_min||node_coordinate(i,1)==x_min||node_coordinate(i,2)==y_max||node_coordinate(i,1)==x_max)
        node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=function_u(node_coordinate(i,1),node_coordinate(i,2));
%         node_coordinate(i,5)=function_v(node_coordinate(i,1),node_coordinate(i,2));
%         node_coordinate(i,6)=function_p(node_coordinate(i,1),node_coordinate(i,2));
[g1,g2]=boundarygt(t,node_coordinate(i,1),node_coordinate(i,2));
 node_coordinate(i,4)=g1;
        node_coordinate(i,5)=g2;
        node_coordinate(i,6)=function_p(node_coordinate(i,1),node_coordinate(i,2),t);
    else
        node_coordinate(i,3)=0;  % 1 = boundary 0= interior
        node_coordinate(i,4)=0;
        node_coordinate(i,5)=0;
        node_coordinate(i,6)=0;
    end
end