function GE=Gauss(A,n)

l=[];
R=[];


for i=1:n-2
    for j=i+1:n
        l(i,j)=A(j,i)/A(i,i)
        R(:,j)=R(:,j)-l(i,j)*R(i,:)
    end
end

%I am confused, but I'll go to office hours on Tuesday; deal?