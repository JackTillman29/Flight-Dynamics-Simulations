function [partialPulse,partialTime] = PartialPulseGrab(fullPulse,fullTime,centerFraction,widthFraction)
    % middle pulse grab
    centerIdx       = round(length(fullTime) * centerFraction );
    grabLength      = round(length(fullTime) * widthFraction);
    partialPulse    = fullPulse((centerIdx-floor(grabLength/2)):(centerIdx+floor(grabLength/2)));
    partialTime     = fullTime( (centerIdx-floor(grabLength/2)):(centerIdx+floor(grabLength/2)));
end