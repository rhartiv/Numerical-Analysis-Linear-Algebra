function [y x]=fdBVPNewton(f,fy,fz,a,b,ya,yb,n)

h=(b-a)/n;

x = linspace(a,b,n+1);
y = linspace(ya,yb,n+1);

J= zeros(n+1,n+1);
J(1,1)=1;
J(n+1,n+1)=1;
F(1)=0;
F(n+1)=0;

tol=10^(-6);
maxdiff=tol+1;
while maxdiff>tol
    for i= 2:n
        ypi=(y(i+1)-y(i-1))/(2*h);
        Pi=.5*h*fz(x(i),y(i),ypi);
        J(i,i-1)=-(1+Pi);
        J(i,i)=2+(h^2)*fy(x(i),y(i),ypi);
        J(i,i+1)=-(1-Pi);
        F(i)=2*y(i)-y(i-1)-y(i+1)+(h^2)*f(x(i),y(i),ypi);
    end
    v=J\F';
    maxdiff=max(abs(v));
    y=y-v;
end