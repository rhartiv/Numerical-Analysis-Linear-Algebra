n=100
a=1; 
b=2;
h=(b-a)/n;
x1=a:h:b;
p=x1(1:n-1).^2;
q=x1(1:n-1);
r=x1(1:n-1).^2;
alpha=0
beta=1
A=zeros(n-1,n-1);

for i=2:n-2
    A(i,i+1)=1-.5*h*p(i);
    A(i,i)=-(2+h^2*q(i));
    A(i,i-1)=1+.5*h*p(i);
end

A(1,1)=-(2+h^2*q(1));
A(1,2)=1-.5*h*p(1);
A(n-1,n-2)=1+.5*h*p(n-1);
A(n-1,n-1)=-(2+h^2*q(n-1));

b=h^2.*r;
b(1)=b(1)-(1+.5*h*p(1))*alpha;
b(end)=b(end)-b(1)-(1+.5*h*p(n-1))*beta;
y1=b/A;
y1=[alpha,y1,beta];

a=1; 
b=2;
%first shoot
s(1)=(beta-alpha)/(b-a);
[x2,y2]=ode45('funsys',[a b],[alpha;s(1);0;1]);

%second shoot
l=length(y2(:,1));
s(2)=s(1)-(y2(l,1)-beta)/y2(l,3);
clear x2;
clear y2;
[x2,y2]=ode45('funsys',[a b],[alpha;s(2);0;1]);

%third shoot
l=length(y2(:,1));
s(3)=s(2)-(y2(l,1)-beta)/y2(l,3);
clear x2;
clear y2;
[x2,y2]=ode45('funsys',[a b],[alpha;s(3);0;1]);
plot(x1,y1,x2,y2(:,1))
