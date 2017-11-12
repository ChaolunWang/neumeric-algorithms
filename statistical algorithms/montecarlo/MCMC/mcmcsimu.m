function [ S ] = mcmcsimu( A, N, times,steps )
    S=zeros(N,N);
    eachelement=zeros(N);
    for k=1:times
        i=randi(N);
        eachelement(i)=eachelement(i)+1;
        q=accept(A(i,:)');
        for j=2:steps
            q=accept(A(q,:)');
        end
        S(i,q)=S(i,q)+1;
    end
    for i=1:N
        S(i,:)=S(i,:)./eachelement(i);
    end
end

