function offset = GenerateTimingJitter(N,samples)
    r = samples*(2*rand(1,N)-1);
    offset = floor(r);
end