function out = vec2idx(in)
    out = zeros(length(in),2);
    out(1,:) = [1 1-1+in(1)];
    for k = 2 : length(in)
        out(k,1) = out(k-1,2)+1;
        out(k,2) = out(k,1) -1 + in(k);
    end
end