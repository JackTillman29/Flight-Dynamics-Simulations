function out = uintArray2Int(a)
    aInternal = uint64(a);
    L = length(a);
    out = uint64(0);
    for k = 1:L
        out = uint64(double(out) + double(2^(k-1)) .* double(aInternal(L-k+1)));
    end
end