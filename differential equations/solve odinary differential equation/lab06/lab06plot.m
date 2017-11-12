% this is the matlab code for plotting the result of lab06 

display('Question1_1 :');               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot FE method result

data = fopen('q11result.txt','r');                              %read data from q11result.txt file
A = fscanf(data ,'%f');
fclose(data);

M=161;                                                          %initialization
T =zeros(M);
X = zeros(M);
Y = zeros(M);
E = zeros(M);
Err = zeros(M);

for i= 1:M                                                       %separate data into array X, T, E, Err and Y
    T(i)=A((i-1)*5+1);
    X(i)=A((i-1)*5+2);
    Y(i)=A((i-1)*5+3);
    E(i)=A((i-1)*5+4);
    Err(i)=A((i-1)*5+5);
end

figure;
plot(T, X, 'r', T, Y, 'b');                                      %plot the data

figure;
plot(X, Y, 'r');

figure;
plot(T, E, 'r');

figure;
plot(T, Err, 'r');   

display('Question1_2 :');                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  plot BE mehtod result

data = fopen('q12result.txt','r');                                 %read data from q12result.txt file
A = fscanf(data ,'%f');
fclose(data);


for i= 1:M                                                         %separate data into array X, T, E, Err and Y
    T(i)=A((i-1)*5+1);
    X(i)=A((i-1)*5+2);
    Y(i)=A((i-1)*5+3);
    E(i)=A((i-1)*5+4);
    Err(i)=A((i-1)*5+5);
end

figure;
plot(T, X, 'r', T, Y, 'b');                                         %plot the data

figure;
plot(X, Y, 'r');

figure;
plot(T, E, 'r');

figure;
plot(T, Err, 'r');

display('Question1_3 :');                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot IT method

data = fopen('q13result.txt','r');                                   %read data from q13result.txt file
A = fscanf(data ,'%f');
fclose(data);


for i= 1:M                                                           %separate data into array X, T, E, Err and Y
    T(i)=A((i-1)*5+1);
    X(i)=A((i-1)*5+2);
    Y(i)=A((i-1)*5+3);
    E(i)=A((i-1)*5+4);
    Err(i)=A((i-1)*5+5);
end

figure;
plot(T, X, 'r', T, Y, 'b');                                            %plot the data

figure;
plot(X, Y, 'r');

figure;
plot(T, E, 'r');
axis([0 35 0 1]);

figure;
plot(T, Err, 'r');

display('Question1_4 :');                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot LF method

data = fopen('q14result.txt','r');                                      %read data from q14result.txt file
A = fscanf(data ,'%f');
fclose(data);


for i= 1:M                                                              %separate data into array X, T, E, Err and Y
    T(i)=A((i-1)*5+1);
    X(i)=A((i-1)*5+2);
    Y(i)=A((i-1)*5+3);
    E(i)=A((i-1)*5+4);
    Err(i)=A((i-1)*5+5);
end

figure;
plot(T, X, 'r', T, Y, 'b');                                             %plot the data

figure;
plot(X, Y, 'r');

figure;
plot(T, E, 'r');

figure;
plot(T, Err, 'r');



