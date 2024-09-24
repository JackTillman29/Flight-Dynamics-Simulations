function y = WaveformRepeatForSamples(u,N)
    % create sized output
    y = zeros(1,N,class(u));
    
    % how many full copies fit??
    nFull = floor(N/length(u));
    
    % how many remaining samples??
    nRem = N-nFull*length(u);
    
    iStart = 1;
    
    % place full copies
    for k = 1 : nFull
        iStop = iStart + length(u) - 1;
        y(iStart:iStop) = u;
        iStart = iStop + 1;    
    end
    
    % place remainder
    y(iStart:end) = u(1:nRem);
    
    

end