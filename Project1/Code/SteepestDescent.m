function [x,k]=SteepestDescent(A,b)

%Find how big it is:
l=length(A);
%Set initial guess
x0=zeros(l,1);
%Set tolerance
tol=10^-6;
%Iitialize iteration count
k=1;
%Loop/Iterate
while norm(b-A*x0) > tol
    %Eliminating large iterations
    if k>=2000
        x='Imaginary array';
        k='Greater than 2000';
        break
    end
    %Setting the v
    v=b-A*x0;
    %Finding the alpha value
    a=v'*v/(v'*A*v);
    %Setting iterative x values
    x=x0+a*v;
    %Revaluing initial
    x0=x;
    %Counting iterations
    k=k+1;
end