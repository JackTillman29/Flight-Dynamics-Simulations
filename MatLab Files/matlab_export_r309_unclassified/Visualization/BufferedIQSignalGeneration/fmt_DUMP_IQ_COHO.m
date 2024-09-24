function out = fmt_DUMP_IQ_COHO(in)
% Data format: double(i,q,i,q...)
endIQ = length(in)/2;
    out = [ in(1:2:endIQ)' + 1i*in(2:2:endIQ)'  ...
            in(endIQ+1:2:end)' + 1i*in(endIQ+2:2:end)' ];
    %out_coho = [ in(endIQ+1:2:end)' + 1i*in(endIQ+2:2:end)' ];
end