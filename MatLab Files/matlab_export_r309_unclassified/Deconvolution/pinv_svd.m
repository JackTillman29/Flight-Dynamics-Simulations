function [Ainv,U,S,V] = pinv_svd(A,N)

if(~exist('N'))
    N = min(size(A));
end

% SVD approach
[U,S,V]=svd(A);

Si = 0*S;
for k = 1 : N
    Si(k,k) = 1.0 / S(k,k);
end


Ainv = (U*Si*V')';
end