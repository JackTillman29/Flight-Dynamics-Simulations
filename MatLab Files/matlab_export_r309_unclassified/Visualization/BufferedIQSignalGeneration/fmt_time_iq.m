function out = fmt_time_iq(in)
% Data format: double(time),double(i,q,i,q...)
    rows = length(in)/3;
    out = [ ...
        in(1:rows)' in((rows+1):2:end)' + 1i*in((rows+2):2:end)' ];
end