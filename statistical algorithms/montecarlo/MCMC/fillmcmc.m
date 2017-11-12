function [ A ] = fillmcmc( In, n )
    A=zeros(n,n);
    for i=1:size(In)
        A(In(i,1),In(i,2))=In(i,3);
    end

end

