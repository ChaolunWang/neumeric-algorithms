clear all; close all; clc

n=2; %number of nodes
N=100;%number of updations
times=1000;
In=[1 2 0.2; 2 1 0.7;1 1 0.8;2 2 0.3];

A=fillmcmc(In, n);

S=A^N;
S
%[V, D]=eig(A);
%V
%D
%[U,T] = schur(A);
%T=T^N;
%for i=1:n
%    T(i,i)=T(i,i)^N;
%end

%Stable=U*T*U';
%Stable
 simulate=mcmcsimu( A, n, times, N );
 simulate
