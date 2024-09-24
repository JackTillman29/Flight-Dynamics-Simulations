function [u,varargout] = gen_pd_impulse_multi(u,rng,prf,Ts,amp,LO,RF,doppler)
% u is pre-allocated array of dwell samples
% rng is the true range of the signal return (m)
% prf (Hz) - used for ambiguous range calculation
% Ts - sample time (s)
% amp - impulse amplitude (v) (drives pulse amplitude)
% LO - Hz, applies appropriate phase shift to impulse terms such
%      that after convolution, signal appears as interrupted CW
%      (e.g. coherent)
% doppler - Hz 
% Optional outputs:
%    varargout(1) = first target sample index into u
%    varargout(2) = timing error of first sample (if TOF does not line up
%                   perfectly with Ts)
%
% NOTE: This function can be called multiple times for multiple 
%       targets as it **adds** to the existing impulse terms. Be
%       sure to complete the impulse vector creation before
%       convolving with the transmitted waveform!!
%
% CAVEATS:
%   1.) Doppler shifts are not captured over the pulse duration; only
%       a constant phase shift. It is not possible (or is it??) to capture
%       a dynamic Doppler over the pulse duration using this impulse method
%       by the nature of the method. This caveat is important when using
%       pulse compression pulses, since Doppler shifts over the pulse
%       duration impact how the pulse is processed by a matched filter.
%           ** possible solution would be design a filter that is convolved
%              with the impulse function to produce the desired doppler
%              shift after the filtered impulse function is convolved with
%              the pulse kernel
%   2.) (NO LONGER A CAVEAT, LEFT HERE SO THAT OTHERS CAN LEARN, SOLUTION
%        DETAILED BELOW HAS BEEN IMPLEMENTED)
%       In general, this method will not produce a pulse train that is
%       coherent with a free-running coherent oscillator.
%       What is meant by "coherent"?
%       If a CW tone, produced with the same starting phase of the pulse
%       kernel that is convolved with the impulse function , is compared to
%       the resultant pulse train derived from the impulse function from
%       this routine, there will be phase discontinuities observed in the
%       impulse-convolution pulse train.
%           ** possible solution is to compute the phase error here
%              internally for each pulse and apply a correction to the 
%              phase of each impulse.
% 
% 
% K. Sawmiller | J. Hole 2019

% compute the equivalent # of samples associated with a PRI (or unambiguous
% range)
arng_samples = round((1/prf) / Ts);

% determine the first sample where a return will appear due to
% time-of-flight
first_rng_sample = max(floor(rng / 150e6 / Ts),1);  % using floor so it's always "early"
range_error = Ts * first_rng_sample * 150e6 - rng;
phase_error = 2*pi*range_error / (3e8 / RF);


numIn = numel(rng);
for k = 1:numIn
    
%     if(mod(k,numIn/20)==0)
%         clc
%         disp(['  building impulse fcn: ',num2str(k/numIn*100),' %'])
%     end
    
    % compute the time vector associated with the targeted impulse samples
    tTimeVec = ((first_rng_sample(k):arng_samples:length(u))-1) * Ts;

    % compute the phase-offset vector associated with the targeted impulse
    % samples. This is needed because the reference transmit waveform starting
    % phase will repeat at each impulse point, and must be "rotated" to address
    % the true phaseing of a real signal
    tPhaseVec = exp(1i*(2*pi*(LO+doppler(k))*tTimeVec));

    % apply sampling phase error
    tPhaseVec = exp(1i*phase_error(k)) .* tPhaseVec;

    % Sum these impulse terms with whatever was there prior
    u(first_rng_sample(k):arng_samples:end) = ...
        u(first_rng_sample(k):arng_samples:end) + amp(k) .* ...
          tPhaseVec;

end

if(nargout > 1)
    varargout{1} = first_rng_sample;
end

end