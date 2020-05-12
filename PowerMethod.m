%Power Method
function [k,u]=PowerMethod(A,x)

n=length(x);
tol=10^-4;
N=25;
k=1;
xp=min(find(abs(x)==(norm(x,inf))));
x=x/x(xp);
while k<=N
    y=A*x;
    yp=min(find(abs(y)==(norm(y,inf))));
    u=abs(y(yp));
    if yp==0
        break
    end
    error=norm((x-(y/u)),inf);
    x=y/u;
    if error<tol
        break
    end
    k=k+1;
end

    
