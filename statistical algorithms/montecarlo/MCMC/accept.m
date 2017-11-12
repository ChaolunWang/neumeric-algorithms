function [ n ] = accept( b)%given the input vector as the possibilities,out put a chosen one
    accumulate=b;
    for i=1:size(b);
        if i>1
            accumulate(i)=accumulate(i)+accumulate(i-1);
        end
    end
    
    q=rand(1);
    for i=1:size(b);
        if accumulate(i)>q
            n=i;
            break;
        end
    end 

end

