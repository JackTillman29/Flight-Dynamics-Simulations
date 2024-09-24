function out = fmt_DUMP_IQ_ONLY(in)
% Data format: double(i,q,i,q...)
    out = [ in(1:2:end)' + 1i*in(2:2:end)' ];
end