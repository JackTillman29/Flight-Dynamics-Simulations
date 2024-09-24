function out = fmt_DUMP_RDM(in)
% Data format: double(i,q,i,q...)
    idx = [1 15 2 16 3 17 4 18 5 19 6 20 7 21 8 22 9 23 10 24 11 25 12 26 13 27 14];
    out1 = reshape(in,27,18);
    out = out1(idx,:);
end