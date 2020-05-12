% Nonlinear

f=@(x,y,z)exp(y)-(1/x)*z
fy=@(x,y,z)exp(y)
fz=@(x,y,z)-(1/x)
a=0.1
b=0.125
ya=-1
yb=1
n=100

[y1 x1]=fdBVPNewton(f,fy,fz,a,b,ya,yb,n)


%first shoot
s(1)=(yb-ya)/(b-a);
[x2,y2]=ode45('funsys',[a b],[ya;s(1);0;1]);

%second shoot
l=length(y2(:,1));
s(2)=s(1)-(y2(l,1)-yb)/y2(l,3);
clear x2;
clear y2;
[x2,y2]=ode45('funsys',[a b],[ya;s(2);0;1]);

%third shoot
l=length(y2(:,1));
s(3)=s(2)-(y2(l,1)-yb)/y2(l,3);
clear x2;
clear y2;
[x2,y2]=ode45('funsys',[a b],[ya;s(3);0;1]);
plot(x2,y2(:,1),x1,y1(60,:))
