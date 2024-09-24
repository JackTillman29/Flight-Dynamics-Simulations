function [Ainv,U,S,V] = pinv_svd_ext(U,S,V,N)

if(~exist('N'))
    N = min(size(U));
end

% SVD approach
%[U,S,V]=svd(A); % now externally passed in

Si = 0*S;
for k = 1 : N
    Si(k,k) = 1.0 / S(k,k);
end


Ainv = (U*Si*V')';
end