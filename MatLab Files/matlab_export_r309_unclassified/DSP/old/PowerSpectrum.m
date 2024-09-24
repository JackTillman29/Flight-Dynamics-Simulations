function ps = PowerSpectrum(frq,cpsd,fbins)
    % frq  - frequency vector of cumulative power spectrum
    % cpsd - Cumulative power over frequency vector frq
    % fbins - Discrete bins of interest
    %   form: [0 1; 1 2; 2 3; 3 4];
    %         Bin 1: All power in 0Hz to 1Hz region
    %         Bin 2: All power in 1Hz to 2Hz region
    %         Bin 3: All power in 2Hz to 3Hz region
    %         Bin 4: All power in 3Hz to 4Hz region

    nBins = size(fbins,1);
    
    ps = zeros(nBins,1);
    
    for k = 1 : nBins
        lowval = find(frq >= fbins(k,1),1,'first');
        hival  = find(frq <= fbins(k,2),1,'last');
        lowpow = cpsd(lowval)
        hipow  = cpsd(hival)
        ps(k) = cpsd(hival)-cpsd(lowval);
    end
    
end

