function y = shiftdata(u,n,option)
% function y = shiftdata(u,n,option)
% u = vector data
% n = # of elements to shift (+ or -)
% option = 'zeropad','circular','expand'
n = round(n);

need_transpose = 0;
if(size(u,1) > size(u,2))
    need_transpose = 1;
    u = u.';
end

    if(abs(n) > length(u))
        error('Shift requested is larger than input!');
    end
        

    switch(lower(option))
        case 'zeropad'
            if(n >= 0)
                y = [zeros(1,n) u(1:(end-n))];
            else
                y = [u((-n+1):end) zeros(1,-n)];
            end
        case 'circular'
            if(n >= 0)
                y = [u((end-n+1):end) u(1:(end-n))];
            else
                y = [u((-n+1):end) u(1:(-n))];
            end
        case 'expand'
            if(n >= 0)
                y = [zeros(1,n) u];
            else
                y = [u zeros(1,n)];
            end
        otherwise
            error(['Bad option: ' option ...
                '. Must be either zeropad or circular']);
    end
    
    if(need_transpose)
        y = y.';
    end
    
end