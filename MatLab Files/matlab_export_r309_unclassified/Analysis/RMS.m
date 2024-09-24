% Updated to handle both real and complex inputs properly
function rms_value = RMS(input)
    rms_value = sqrt(mean(input.*conj(input)));
    if(~(isreal(input)))
        rms_value = rms_value * sqrt(2) / 2;
    end
end