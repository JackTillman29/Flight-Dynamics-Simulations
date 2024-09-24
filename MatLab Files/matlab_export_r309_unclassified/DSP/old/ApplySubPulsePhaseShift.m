function varargout = ApplySubPulsePhaseShift(varargin)
    pulse = varargin{1};
    phaseShifts = varargin{2}; %real phase shift keys
    
	% determine the number of sub-pulses
    nLengthSubPulses = floor(length(pulse)./length(phaseShifts));
    nSubPulses       = floor(length(pulse)./nLengthSubPulses);
    
    for i = 1:nSubPulses
        starti = (i-1)*nLengthSubPulses+1;
        stopi  = starti + nLengthSubPulses-1;
        complexPhaseShifter = exp(1j*phaseShifts(i));
        pulse(starti:stopi) = complexPhaseShifter*pulse(starti:stopi);
    end
    if ( stopi < length(pulse) )
        pulse = pulse(1:stopi);
        warning(['Truncated input pulse because the subpulses did not fit evenly']);
    end
    varargout{1} = pulse;
end