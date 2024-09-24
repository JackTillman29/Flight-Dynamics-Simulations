function w = window_cheby2(N,AdB)
    w = zeros(1,N);
    r_inv = 10^(-AdB/20);
    N_inv = 1/N;
    x0 = cosh( 1/(N-1) * acosh(r_inv) );
    
    for n = 1 : N
        sum_term = 0;
        for ki = 1 : floor((N-1)/2)
            T_term = (x0 * cos(ki*pi/N));
            sum_term = sum_term + ...
                T(N-1,T_term)*cos(2*pi*n*ki/N);
        end
        
        w(n) = N_inv * ( ...
            r_inv + ...
            2*sum_term ...
            );
    end
end

function chebpoly = T(n,x)
    if(abs(x) <= 1)
        chebpoly = cos(n*acos(x));
    else
        chebpoly = cosh(n*acosh(x));
    end
end