function [num,den] = tfmult(num1,den1,num2,den2)

num = polymult(num1,num2);
den = polymult(den1,den2);
end