classdef LFMPulse < Pulse
    properties 
        sweepRate; % hz/sec
    end
    methods
        function p = LFMPulse(sweepRate,varargin)
            p = p@Pulse(varargin{:});
            if ( nargin == 0 )
               p.sweepRate = 0;
            elseif ( nargin == 1 ) 
               p.sweepRate = sweepRate;
            else
               p.sweepRate = sweepRate;
            end
            
            % get max frequency for sampling
            fmax = max(p.carrier,p.carrier + p.pulsewidth * p.sweepRate);
            p.Fs = p.oversample * fmax;
            p.Ts = 1.0 / p.Fs;
            p.generateTDSignal;
        end
        
        function generateTDSignal(p)
            p.timeVector      = 0:p.Ts:p.pulsewidth;
            freqLFM = interp1([0 p.pulsewidth],[p.carrier p.carrier + p.pulsewidth * p.sweepRate],p.timeVector);
            p.amplitudeVector = p.windowfcn(length(p.timeVector))' .* p.peakAmplitude .* ....
                exp(1j*2*pi*freqLFM.*p.timeVector);
        end
        function pf = frequencyOffsetPulse(p,df)
            pf = LFMPulse(0);
            p.copyProperties(pf);
            pf.carrier = p.carrier + df;
            pf.Fs = p.Fs;
            pf.Ts = p.Ts;
            pf.generateTDSignal;
        end
    end
end