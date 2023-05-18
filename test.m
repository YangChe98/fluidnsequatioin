A=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            A(i,j,k)=i*100+j*10+k;
        end
    end
end

B=[1 2 3 ].';
B1=repmat(B,1,3,3);
sumA=sum(B1.*A,1);
sumA1=zeros(3,3);
for k=1:3
    sumA1=sumA1+B(k)*A(:,:,k);
end
