function [quantized,bitsNeeded] = ADC(inputSignal,nBits,maxV)
    % 
    if(~exist('maxV','var'))
        vRms = sqrt(mean(inputSignal.^2));
        maxV = 1.25 * 2.0 * vRms / sqrt(2); % Vpeak = 2Vrms / sqrt(2), 25% padding
        disp(['Auto calculating Vmax: ' num2str(maxV)]);
    end

    % the number of discrete values available is 2^N, where N = # of bits
    nValues = 2.0 ^ nBits;
    
    % this model assumes the ADC can sample both (+/-) voltages
    vPerBit = maxV / ((nValues-1)/2);
    
    % build map of voltage to bits
    bitMap = ((0:(nValues-1)) * vPerBit) - maxV;
    
    % temp
    temp = single(0.0 * inputSignal);
    
    for k = 1 : length(bitMap)
        idx = find(inputSignal >= (bitMap(k)-0.5*vPerBit));
        temp(idx) = k-1;
    end
    
    bitsNeeded = length(inputSignal) * nBits;
    %disp(['Storage requirement: ' num2str(length(inputSignal) * nBits / 1024) 'Kb']);
    
    % map signal back to voltage
    quantized = temp * vPerBit - maxV;
    
    
end