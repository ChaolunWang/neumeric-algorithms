clear all; close all; clc

n=11;
dt=0.05;
q=1;
D=4;

%calculate h
h=1/n;
xspan=0:h:1;

%build S
e1=ones(n-1,1);
S=spdiags([-1*e1 2*e1 -1*e1], [-1 0 1], n-1, n-1);
S=(1/h).*S;

%build M
M=spdiags([e1 4*e1 e1], [-1 0 1], n-1, n-1);
M=(h/6).*M;

%fxt
fxt= @(x,t)((1-D*pi^2).*(exp(t)).*sin(pi.*x));
f= @fem1dbase;
gx=@(x,xi,h,t) fxt(x,t).*f(xi,x,h);
bx=@(xi,h,t) (gx(xi-1/2*h,xi,h,t)+gx(xi+1/2*h,xi,h,t)).*h;

%t=0;
%b=bx(xspan(1,2:n+1), h, t)';
%initial condition
ux0=@(x) sin(pi.*x);
gx0=@(x,xi,h) ux0(x).*f(xi,x,h);
bx0=@(xi,h) (gx0(xi-1/2.*h,xi,h)+gx0(xi+1/2.*h,xi,h)).*h;
b0=bx0(xspan(1,2:n), h)';
eta=M\b0;


[L,U]=lu(M-D.*S.*dt);
for ti=dt:dt:q
    b=bx(xspan(1,2:n), h, ti)';
    eta=U\(L\(b.*dt+M*eta));
end

y=zeros(n+1,1);
eta



for i=1:n-1
    y(i+1)=0;
    for j=1:n-1
        y(i+1)=y(i+1)+eta(j)*fem1dbase(j*h,i*h,h);
    end

end

uxt=@(x) exp(q)*sin(pi.*x);

plot(xspan, y' ,'.');
hold;
plot(xspan, uxt(xspan), 'g');