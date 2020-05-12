% Inverse Power Method
function [k,y]=InvPowerMethod(A,x)

n=length(x);
tol=10^-6;
N=3;
k=1;
q=(x.')*A*x/((x.')*x)
xp=min(find(abs(x)==(norm(x,inf))));
x=x/abs(x(xp));
while k<=N
    y=inv(A-q*eye(n))*x;
    yp=min(find(abs(y)==(norm(y,inf))));
    u=y(yp);
    if yp==0
        break
    end
    error=norm((x-(y/u)),inf);
    x=y/y(yp);
    if error<tol
        u=(1/u)+q
        break
    end
    k=k+1;
end

    
