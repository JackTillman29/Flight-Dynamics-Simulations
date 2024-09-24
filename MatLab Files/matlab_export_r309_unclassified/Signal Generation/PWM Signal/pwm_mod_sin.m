function [outScalePri,outScalePw] = pwm_mod_sin(inFrac,ampl,cycles,ampl_pw,cycles_pw)
    % inFrac expects to go from 0 to 1 over the waveform generation
    % think of it as "fraction of technique period"
    outScalePri = ampl*sin(2*pi*cycles*inFrac)+1;
    outScalePw = ampl_pw*sin(2*pi*cycles_pw*inFrac)+1;
    
end