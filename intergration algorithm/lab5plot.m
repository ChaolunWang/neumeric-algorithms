% this is the matlab code for plotting the result of lab5

display('Question1 :');

data = fopen('q12result.txt','r');                     %read data from q12result.txt file
A = fscanf(data ,'%f');
fclose(data);

X = zeros(101);
Y = zeros(101);
n = zeros(101); 
X2 = zeros(5);
Y2 = zeros(5);
X3 = zeros(11);
Y3 = zeros(11);
for i= 1:101                                           %separate data into array X, n and Y
    X(i)=A((i-1)*3+1);
    Y(i)=A((i-1)*3+2);
    n(i)=A((i-1)*3+3);
end

figure;
plot(X,Y);                                             %plot the k(x)

display('Question2 :');
figure;
plot(X, n);                                             %plot the value of n for each x

display('Question3 :');

data = fopen('q3_1.txt','r');                         %read data from q3-1.txt file
B1 = fscanf(data ,'%f');
fclose(data);

for i= 1:101                                           %separate data into array X and Y
    X(i)=B1((i-1)*2+1);
    Y(i)=B1((i-1)*2+2);
end

figure;
plot(X,Y);                                             %plot the intergrand

data = fopen('q3_2.txt','r');                        %read data from q3-2.txt file
B2 = fscanf(data ,'%f');
fclose(data);

for i= 1:5                                           %separate data into array X and Y
    X2(i)=B2((i-1)*2+1);
    Y2(i)=B2((i-1)*2+2);
end

figure;
plot(X2,Y2, 'p');                                       %plot the absissa and weight



data = fopen('q3_3.txt','r');                         %read data from q3-3.txt file
B3 = fscanf(data ,'%f');
fclose(data);

for i= 1:11                                           %separate data into array X and Y
    X3(i)=B3((i-1)*2+1);
    Y3(i)=B3((i-1)*2+2);
end

figure;
plot(X3 ,Y3, 'p');                                         %plot the absissa and weight




