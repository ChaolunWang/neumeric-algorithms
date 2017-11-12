function  [err]=q4function(I)
%This function can calculate the error of bspline approximation based on different number of bins,
%it will read the data form xy.dat file, use bspline package function(bspline_basismatrix) 
%to get the basis matrix and use the cholesky decomposition to solve the linear least square
% problem.
%parameter:
%           I: the number of bins of approximation
%output:
%           err: the error of approximation
    data = fopen('xy.dat','r');                %read data form xy.dat file
    A = fscanf(data ,'%f');
    fclose(data);

    for i= 0:499                               %separate the data into array X and Y
        X(i+1)=A(i*2+1);
        Y(i+1)=A(i*2+2);
    end
    P=2; M=I+P*2;
    
    for i= 1:M+1                               %build knod vector
        if i<=P
            T(i)=X(1);
        elseif i>M+1-P;
            T(i)=X(500);
        else
            T(i)=(X(500)-X(1))/I*(i-P-1);
        end   
    end
    A= bspline_basismatrix(P+1,T,X);          %call bspline_basismatirx to get basis matrix

    AT=transpose(A);                          %use cholesky decomposition to solve the linear least square problem
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
    
    for i=1:500                               %build the approximation function
        Y0(i)=0;
        for j=0 : numel(x_c)-1
            Y0(i)=Y0(i)+ bspline_basis(j,P+1,T,X(i))*x_c(j+1);
        end
    
    end

    err=0;                                     %calculate error
    for i=1:500   
        err=err+(Y(i)-Y0(i))^2;
    end
end
