
%SOR
function [x,k]=SOR(A,b,omega)
%Find how big it is:
n=length(b);
%Set tolerance
tol=10^-6;
%Set initial guess
x=zeros(n,1);
%Create D,U,L,T, and c
D=diag(A).*eye(length(A));
U=-triu(A)+D;
L=-tril(A)+D;
T=inv(D-omega*L)*(omega*U+(1-omega)*D);
c=omega*inv(D-omega*L);
%Set initial error
error=norm(A*x-b,inf);
%Iitialize iteration count
k=0;
%Loop/Iterate
while error>tol
    %Eliminating large iterations
    if k>=2000
        x='Imaginary array';
        break
    end
    %Setting iterative x values
    x=T*x+c*b;
    %Editing the error
    error=norm(A*x-b,inf);
    %Counting iterations
    k=k+1;
end
end 