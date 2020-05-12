function [x,tol]=GaussSeidelIteration(A,b,k)
n=length(b);
x=zeros(n,1);
D=diag(A).*eye(length(A));
U=-triu(A)+D;
L=-tril(A)+D;
T=inv(D-L)*U;
c=inv(D-L)*b;
error=norm(A*x-b,inf);
for i=1:k
    x=T*x+c;
    tolerance=10^floor(log10((max(abs(A*x-b)))));
end
tol=tolerance;
end 
