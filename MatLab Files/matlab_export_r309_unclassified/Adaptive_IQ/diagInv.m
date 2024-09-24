function Uinv = diagInv(U,N)
% This function takes the diagonal elements, inverts them element-wise
% and then returns them as a new diagonal matrix
Uinv = diag(1.0 ./ diag(U));

if(exist('N'))
    
    for k = (N+1):length(U)
        Uinv(k,k) = 0;
    end
end
end