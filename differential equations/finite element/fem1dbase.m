function [ v ] = fem1dbase( xi, x, h)
%1D base function

    if (xi-h<=x) & (x<xi)
        v=(x-(xi-h))/h;
    elseif (xi<=x) & (x<xi+h)
        v=(xi+h-x)/h;
    else
        v=0;
    end 
end

