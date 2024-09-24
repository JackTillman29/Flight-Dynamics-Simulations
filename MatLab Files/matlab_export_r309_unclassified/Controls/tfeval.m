function output = tfeval(num,den,k,input)
    output = polyval(k*num,input) ./ polyval(den,input);
end