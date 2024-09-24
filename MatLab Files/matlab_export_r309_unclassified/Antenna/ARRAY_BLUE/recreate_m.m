function Mout = recreate_m(M)

    [V,D]=eig(M);

    N = length(diag(D));

    Mout = zeros(N);
    
    for k = 1 : N
        Mout = Mout + ...
            D(k,k) * V(:,k) * V(:,k)';
    end

end