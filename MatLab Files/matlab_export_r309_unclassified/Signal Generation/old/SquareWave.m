function [output,ampitude_terms,frequency_terms] = SquareWave(frequencyHz,timeVector,nTerms,dcOffset,atten)
    output = zeros(1,length(timeVector));

    amplitude_terms     = zeros(1,nTerms+1);
    frequency_terms     = zeros(1,nTerms+1);
    amplitude_terms(1)  = dcOffset;
    frequency_terms(1)  = 0.0;
    output              = output + dcOffset;
    
    for k = 1:nTerms
        amplitude_terms(k+1) = atten * (4.0 / pi * 1.0 / (2.0*k - 1.0));
        frequency_terms(k+1) = (2*k-1)*frequencyHz;
        output = output + amplitude_terms(k+1) .* sin(2*pi*frequency_terms(k+1)*timeVector);
    end
end