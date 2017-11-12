close all; clear all; clc

TOL=0.000001;
n=200;
dx=20/n;
x=zeros(n,1);
for i=1:n;
   x(i)=i*dx; 
end
b0t=5; b1t=2; b2t=3; b3t=1;
fun= @(t) b0t+b1t.*exp(-1.*(t./b2t).^b3t);

y=fun(x)+ 0.5*randn(size(x));

plot(x,y);
hold on;

B=ones(4,1);
B(2)=4;
B(3)=2;
r=zeros(n,1);
J=zeros(n,4);
n1=0;
n2=-1;
while abs(n1-n2)>TOL
    n2=n1;
    for i=1:n
       r(i)=y(i)-B(1)-B(2)*exp(-1*(x(i)/B(3))^B(4)); 
       J(i,1)=-1;
       J(i,2)=-1*exp(-1*(x(i)/B(3))^B(4));
       J(i,3)=-1*B(2)*B(4)*exp(-1*(x(i)/B(3))^B(4))*x(i)^B(4)*B(3)^(-1*B(4)-1);
       J(i,4)=B(2)*exp(-1*(x(i)/B(3))^B(4))*(x(i)/B(3))^B(4)*log(x(i)/B(3));
    end
    n1=norm(r);
    if cond(J'*J)>10000
        break;
    end
    B=B-inv(J'*J)*J'*r;
end
B
cond(J'*J)
fun= @(t) B(1)+B(2).*exp(-1.*(t./B(3)).^B(4));
fplot(fun, [0 20])
hold off;