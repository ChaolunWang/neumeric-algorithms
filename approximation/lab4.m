%This is the main code for lab4, in which the bspline method was used for
%approximate the noizy data. The bspline package was used for the bspline
%algorithm. Be sure to inclue the bspline folder into the path before running
%the code.

%programmed by: Chaolun Wang
%           at: 03/05/2016


%code for question 1
display('Question1 :');

data = fopen('xy.dat','r');                            %read data from xy.data file
A = fscanf(data ,'%f');
fclose(data);

for i= 0:499                                           %separate data into array X and Y
    X(i+1)=A(i*2+1);
    Y(i+1)=A(i*2+2);
end

figure;
plot(X,Y);                                             %plot the noizy data
hold on;


P=2; M=14; I=10;

%create knot vector T:
for i= 1:M+1
    if i<=P
        T(i)=X(1);
    elseif i>M+1-P;
        T(i)=X(500);
    else
        T(i)=(X(500)-X(1))/I*(i-P-1);
    end   
end


A= bspline_basismatrix(P+1,T,X);                      %using bspline_basismatrix function to get basis matrix

AT=transpose(A);
x_A=linsolve(AT*A,AT*transpose(Y));                   %solve linear system to find the control points


for j=0 : numel(x_A)-1                                %plot basis functions
    for i=1:500
        Y0(i)=bspline_basis(j,P+1,T,X(i))*x_A(j+1);
    end
    if mod(j,3)==0
        plot(X, Y0, 'g');
    elseif mod(j,3)==1
        plot(X, Y0, 'r');
    else
        plot(X, Y0, 'y');
    end
end

for i=1:500                                            %combine basis functions and plot the approximation function
    Y0(i)=0;
    for j=0 : numel(x_A)-1
        Y0(i)=Y0(i)+ bspline_basis(j,P+1,T,X(i))*x_A(j+1);
    end
    
end

plot(X,Y0,'k');                                    

hold off;


%code for question 2
display('Question2 :');
conditionNumber=cond(A);                                    %calculate condition number
conditionNumber

%Cholesky decomposition

tic;                                                         %calculate computational time

C=AT*A;
L = chol(C,'lower');
LT = chol(C,'upper');
d=AT*transpose(Y);
opts.LT = true;
z=linsolve(L,d,opts);
opts.LT = false;
opts.UT = true;
x_c=linsolve(LT,z,opts);
opts.UT = false;


wtime=toc;                                                     %calculate computational time
fprintf ( 1, 'Elapsed time = %f\n', wtime );

for i=1:500
    Y0(i)=0;
    for j=0 : numel(x_c)-1
        Y0(i)=Y0(i)+ bspline_basis(j,P+1,T,X(i))*x_c(j+1);
    end
    
end

err=0;                                                            %calculate error
for i=1:500   
    err=err+(Y(i)-Y0(i))^2;
end
fprintf ( 1, 'Least squared error is %f\n', err );

%QR decomposition

tic;                                                             %calculate computational time

[Q,R] = qr(A);
d=transpose(Q)*transpose(Y);
opts.UT = true;
x_q=linsolve(R,d,opts);
opts.UT = false;


wtime=toc;                                                        %calculate computational time
fprintf ( 1, 'Elapsed time = %f\n', wtime );

for i=1:500
    Y0(i)=0;
    for j=0 : numel(x_q)-1
        Y0(i)=Y0(i)+ bspline_basis(j,P+1,T,X(i))*x_q(j+1);
    end
    
end

err=0;                                                             %calculate error
for i=1:500   
    err=err+(Y(i)-Y0(i))^2;
end
fprintf ( 1, 'Least squared error is %f\n', err );


for i=1:5                                                           %print result
    q2(i)=i*10;
    [err1(i),err2(i),time1(i),time2(i)]=q2function(i*10);
    fprintf ( 1, 'When nbin= %d\n', i*10 );
    fprintf ( 1, 'Cholesky decomposition least squared error  is %f\n', err1(i) );
    fprintf ( 1, 'time spent: %f\n', time1(i) );
    fprintf ( 1, 'QR algorithm least squared error  is %f\n', err2(i) );
    fprintf ( 1, 'time spent: %f\n', time2(i) );
end
figure;                                                              %plot error related to nbin
plot(q2,err1,'p');
hold on;
plot(q2,err2,'r');
hold off;

figure;                                                               %plot computational time related to nbin
plot(q2,time1,'g');
hold on;
plot(q2,time2,'r');
hold off;

%code for question 3
display('Question3 :');
for p=0:5
    q3(p+1)=p;
    err3(p+1)=q3function(p);                                          %call q3function to solve question 3
    fprintf ( 1, 'When p= %d', p );
    fprintf ( 1, ', least squared error  is %f\n', err3(p+1) );
end
figure;
plot(q3,err3);


%code for question 4
display('Question4 :');
for i=1:10
    q4(i)=i*10;
    err4(i)=q4function(i*10);                                          %call q4function to solve question 4
    fprintf ( 1, 'When nbin= %d', i*10 );
    fprintf ( 1, ', least squared error  is %f\n', err4(i) );
end
figure;
plot(q4,err4);
