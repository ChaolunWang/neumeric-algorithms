%this is the matlab code for lab3 question 2
display('Question2 a:');
D= input('Please enter the dimention of approximation: ');            %pomp the user for the input
Filename= input('Please enter the file name: ', 's');                 %input the specific filename
A=imread(Filename);                                                   %load file

s = svds(double(A), 150);                                             %do partial svd on A, save singular value, change the 150 into D, if other dimension is needed
figure;
plot(s(1:150), 'b'); 


display('Question2 b:');

figure;
[U, S, V]= svd(double(A));
for i=1:5
    D=4*(2^i);                                                        %set the rank to different value
RanknAprox=double(zeros(512,512));
    for n=1:D
        RanknAprox=RanknAprox+U(: , n)*V(: , n)'*S(n , n);            %calculate different rank of approximation
    end
    subplot(2,3,i), imshow(RanknAprox, [0,255]);
end



subplot(2,3,6), imshow(A);


display('Question2 c:');


r=512;

for p=1:512
    term1=0;
    term2=0;

    for i=(p+1):r
        term1=term1+S(i, i)^2;
    end

    for i=1:r
        term2=term2+S(i, i)^2;
    end

    RelativeError2=sqrt(term1/term2);
    
    if (RelativeError2<0.005)                           %break the loop untill the relative error smaller than 0.005
        break;
    end
end

p

Approx=double(zeros(512,512));                          %display the approximation
for n=1:p
    Approx=Approx+U(: , n)*V(: , n)'*S(n , n);
end
figure;
imshow(Approx, [0,255]);


