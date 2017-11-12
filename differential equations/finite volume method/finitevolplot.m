clear all; close all; clc

n=40;
dx=10/n;
v=2;
C=0.4;

dt=C*dx/v;
%dt=0.005;
a=0;
b=10;
k=(b-a)/dt+1;

nt=(b-a)/dt;


A=importdata('solution0.4');
TV=importdata('tv0.4');

for i=1:nt+1
   tsp(i,1)=a+(i-1)*dt; 
end

figure;
 colormap('hot');   % set colormap
 imagesc([-5, 5], [a, b], A);        % draw image and scale colormap to values range
 colorbar;          % show color scale
 
 
 figure;
 plot(tsp, TV);
 
 %C=0.8
 C=0.8;

dt=C*dx/v;
%dt=0.005;
a=0;
b=10;
k=(b-a)/dt+1;

nt=(b-a)/dt;


A=importdata('solution0.8');
TV=importdata('tv0.8');

for i=1:nt+1
   tsp2(i,1)=a+(i-1)*dt; 
end

figure;
 colormap('hot');   % set colormap
 imagesc([-5, 5], [a, b], A);        % draw image and scale colormap to values range
 colorbar;          % show color scale
 
 
 figure;
 plot(tsp2, TV);
 
 
 %C=1
 C=1;

dt=C*dx/v;
%dt=0.005;
a=0;
b=10;
k=(b-a)/dt+1;

nt=(b-a)/dt;


A=importdata('solution1');
TV=importdata('tv1');

for i=1:nt+1
   tsp3(i,1)=a+(i-1)*dt; 
end

figure;
 colormap('hot');   % set colormap
 imagesc([-5, 5], [a, b], A);        % draw image and scale colormap to values range
 colorbar;          % show color scale
 
 
 figure;
 plot(tsp3, TV);
 
 
 %C=1.2
 C=1.2;

dt=C*dx/v;
%dt=0.005;
a=0;
b=10;
k=(b-a)/dt+1;

nt=(b-a)/dt;


A=importdata('solution1.2');
TV=importdata('tv1.2');

for i=1:nt+1
   tsp4(i,1)=a+(i-1)*dt; 
end

figure;
 colormap('hot');   % set colormap
 imagesc([-5, 5], [a, b], A);        % draw image and scale colormap to values range
 colorbar;          % show color scale
 
 
 figure;
 plot(tsp4, TV);
 
 
 
 %%%%%%%%%%%%%%%%%%%%%part 2 calculate convergence
 
 
 x=[80, 160, 320, 640];
 y=[0.032606, 0.0184734, 0.00988794, 0.00494853];
 
 figure;
 loglog(x,y);
 %the convergent rate for n=160, 320 and 640 are:
 log(y(1,2)/y(1,1))/log(1/2)
 log(y(1,3)/y(1,2))/log(1/2)
 log(y(1,4)/y(1,3))/log(1/2)