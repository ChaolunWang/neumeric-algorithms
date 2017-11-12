clear all; close all; clc

Rn=10000;
rc=1;
n=100;
C=0:2*pi/100:2*pi;
x=rc*cos(C);
y=rc*sin(C);

plot(x,y,'r');
hold;
counter=0;
err=zeros(Rn);

mu=0;%for normal distribution
sigma=0.25;%for normal distribution
%Xn=mu+sigma*rand(1);
for i=1:Rn
    a=-1*rc;
    b=rc;
    
    Xn=a+(b-a)*rand(1);
    Yn=a+(b-a)*rand(1);
    
    if sqrt(Xn^2+Yn^2)<rc
        counter=counter+1;
       plot(Xn,Yn,'.r');
       %hold;
    else
       plot(Xn,Yn,'.b');
       %hold;        
    end
    err(i)=(counter/i*4-pi)/pi;
end

estpi=counter/Rn*4;
figure;
plot(1:1:Rn,err,'b');
estpi
