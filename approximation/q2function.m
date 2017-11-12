function  [err1,err2,time1,time2]=q2function(I)
%this function can get the number of bins(nbin) as input, then get data form xy.dat document
%and output the calculation time and error of Cholesky decomposition and QR decomposition
%parameter:
%           I: number of bins
%output:
%           err1: error of Cholesky method
%           err2: error of QR decomposition
%           time1: time spent by Cholesky method
%           time2: time spent by QR method
    data = fopen('xy.dat','r');                          %read data form file xy.dat
    A = fscanf(data ,'%f');
    fclose(data);

    for i= 0:499                                         %separate the data into array x and y
        X(i+1)=A(i*2+1);
        Y(i+1)=A(i*2+2);
    end

    P=2; M=I+P*2;
    
    for i= 1:M+1                                         %build the knod vector
        if i<=P
            T(i)=X(1);
        elseif i>M+1-P;
            T(i)=X(500);
        else
            T(i)=(X(500)-X(1))/I*(i-P-1);
        end   
    end
    A= bspline_basismatrix(P+1,T,X);

    AT=transpose(A);
    tic;                                                      %calculate time spent
    %Cholesky method
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
    time1=toc;                                              %calculate time spent

    for i=1:500                                             %get the approximation
        Y0(i)=0;
        for j=0 : numel(x_c)-1
            Y0(i)=Y0(i)+ bspline_basis(j,P+1,T,X(i))*x_c(j+1);
        end    
    end

    err1=0;                                                  %calculate error
    for i=1:500   
        err1=err1+(Y(i)-Y0(i))^2;
    end

    tic;                                                     %calculate time spent
    %QR algorithm
    [Q,R] = qr(A);
    d=transpose(Q)*transpose(Y);
    opts.UT = true;
    x_q=linsolve(R,d,opts);
    opts.UT = false;

    time2=toc;                                                %calculate time spent

    for i=1:500
        Y0(i)=0;
        for j=0 : numel(x_q)-1
             Y0(i)=Y0(i)+ bspline_basis(j,P+1,T,X(i))*x_q(j+1);
        end
    end

    err2=0;
    for i=1:500   
        err2=err2+(Y(i)-Y0(i))^2;                            %calculate error
    end
end
