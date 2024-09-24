classdef Pulse < hgsetget
    properties
        pulsewidth;
        windowfcn;
        carrier;
        oversample;
        peakAmplitude;
        Ts;
        Fs;
        timeVector;
        amplitudeVector;
        noisepower;
    end
    methods
        function p = Pulse(varargin)
            
            if (nargin == 1 & isa(varargin{1},'Pulse'))
                varargin{1}.copyProperties(p);
            else
            % set defaults
            p.pulsewidth = 10e-6;
            p.windowfcn  = str2func('@rectwin');
            p.carrier    = 25e6;
            p.oversample = 20;
            p.peakAmplitude = 1;
            p.noisepower = 0;
            for k = 1:2:nargin
                try
                switch varargin{k}
                    case 'pulsewidth'
                        p.pulsewidth = varargin{k+1};
                    case 'amplitude'
                        p.peakAmplitude = varargin{k+1};
                    case 'windowfcn'
                        p.windowfcn  = varargin{k+1};
                    case 'carrier'
                        p.carrier    = varargin{k+1};
                    case 'oversample'
                        p.oversample = varargin{k+1};
                    case 'noisepower'
                        p.noisepower = varargin{k+1};
                    otherwise
                        error(['Invalid property: ' varargin{k}]);
                end
                catch
                    % field is a derived class
                end
                
            end
            p.computeSampleParms;
            
            end
            if ( strcmp(p,'Pulse') == 1 ) % 
                p.generateTDSignal();
            end
        end
        
        function computeSampleParms(p)
            p.Fs = p.oversample * p.carrier;
            p.Ts = 1.0 / p.Fs;
        end
        
        % Copies all member properties to the input object
        function copyProperties(this,new)
            try
                fieldList = fields(this);
                for k = 1:length(fieldList)
                    eval(['new.' fieldList{k} ' = this.' fieldList{k} ';']);
                end
            catch
                %disp('Reverting to base class copyProperties for Pulse');
                fieldList = fields(Pulse);
                for k = 1:length(fieldList)
                    eval(['new.' fieldList{k} ' = this.' fieldList{k} ';']);
                end
            end
                
        end
        
        function generateTDSignal(p)
            p.timeVector      = 0:p.Ts:p.pulsewidth;
            p.amplitudeVector = p.windowfcn(length(p.timeVector))' .* p.peakAmplitude .* ....
                exp(1j*2*pi*p.carrier*p.timeVector) + sqrt(p.noisepower).*exp(1j*2*pi*randn(1,length(p.timeVector)));
            
        end
        
        function [output,kleadlag] = computeAutoCorrelation(p)
            [output,kleadlag] = computeCrossCorrelation(p,p);
        end
        
        function [output,kleadlag] = computeCrossCorrelation(p,s) % s is external sig
            if ( p.Ts ~= s.Ts )
                error('Sample Times Are NOT Consistent! Results will be invalid');
            end
            % Subtract any residual mean (DC) off the signals
            %inputWaveform = inputWaveform - mean(inputWaveform);
            %referenceWaveform = referenceWaveform - mean(referenceWaveform);
            
            % Determine appropriate power of 2
            NFFT = 2^(nextpow2(max([length(p.amplitudeVector) length(s.amplitudeVector)]))+1);
            
            % Compute the FFT of each waveform
            fftIn  = fft(s.amplitudeVector ,NFFT);
            fftRef = fft(p.amplitudeVector ,NFFT);
            
            % Fourier theory says that the FFT of the cross correlation of two signals in the
            % time domain is equal to the product of one signal FFT times the
            % complex conjugate of the other signal FFT
            xcorrFFT = fftIn .* conj(fftRef);
            
            % Determine autocorrelation values for each signal
            inputAutoCorr = ifft(fftIn.*conj(fftIn));
            refAutoCorr   = ifft(fftRef.*conj(fftRef));
            
            % Compute the time domain cross correlation via the IFFT
            xcorrTD = real(ifft(xcorrFFT) ./  ...
                (sqrt(inputAutoCorr(1)) * sqrt(refAutoCorr(1))));
            
            %
            li = length(s.amplitudeVector);
            lr = length(p.amplitudeVector);
            
            xcorrTD = xcorrTD([(end-lr+2):end 1:(end-lr)]);
            
            
            output   = xcorrTD(1:(li+lr-1));
            kleadlag = (0:(li+lr-2)) - round((li+lr-2)/2);
            
            % if no arguments requested, go ahead and plot
            if ( nargout == 0 )
            figure;
            % LeadLag samples are in units of Ts. To convert to (Factor) x
            % (PulseWidth), multiply by Ts/pulsewidth
            plot(kleadlag .* p.Ts ./ p.pulsewidth,[(abs(output'))]);
            % plot(LeadLag ,[(abs(AutoCor'))]);
            axis tight;
            grid on;
            xlabel('Lead/Lag PCT of Pulsewidth');
            ylabel('Correlation');
            title('Pulse Cross Correlation Magnitude');
            end
        end
        
        function plot(this)
            figure;
            
            subplot(2,1,1);
            plot(1e6*this.timeVector, ...
                [real(this.amplitudeVector)' imag(this.amplitudeVector)']);
            axis tight;
            grid on;
            xlabel('Time (\musec)');
            ylabel('Amplitude');
            legend('I','Q','Location','Best');
            title('Pulse Time History');
            
            subplot(2,1,2);
            [AutoCor,LeadLag] = this.computeAutoCorrelation;
            % LeadLag samples are in units of Ts. To convert to (Factor) x
            % (PulseWidth), multiply by Ts/pulsewidth
            plot(LeadLag .* this.Ts ./ this.pulsewidth,[(abs(AutoCor'))]);
            % plot(LeadLag ,[(abs(AutoCor'))]);
            axis tight;
            grid on;
            xlabel('Lead/Lag PCT of Pulsewidth');
            ylabel('Correlation');
            title('Pulse Autocorrelation Magnitude');
        end
        
        function pp = partialPulse(p,fctStart,fctStop)
            iStart = max(round(length(p.timeVector) * fctStart),1);
            iStop  = round(length(p.timeVector) * fctStop);
            pp = Pulse(p);
            pp.timeVector = pp.timeVector(iStart:iStop);
            pp.amplitudeVector = pp.amplitudeVector(iStart:iStop);
        end
        
        function pf = frequencyOffsetPulse(p,df)
            pf = Pulse(p);
            pf.carrier = p.carrier + df;
            pf.generateTDSignal;
        end
        
        function [pri_amp,pri_time] = getInPri(p,dc)
            A = length(p.timeVector);
            n = ( A - dc*A ) ./ dc;
            pri_amp = [p.amplitudeVector zeros(1,n)];
            pri_time = p.Ts * (0:(length(pri_amp)-1));
            
        end
        
        function [A,fv,tv]=ambiguityFunction(p,pwpct,pwpts,frqpct,frqpts)
            % percentages are +/-, so frqpct of 10khz would be +/- 10 khz
            % make sure points are odd length for centering reasons
            if ( mod(pwpts,2) ~= 1 )
                pwpts = pwpts + 1;
            end
            if ( mod(frqpts,2) ~= 1 )
                frqpts = frqpts + 1;
            end

            
            
            df = 2*p.carrier * frqpct / frqpts;
            dt = 2*p.pulsewidth * pwpct / pwpts;
            
            fv  = df * ((1:frqpts) - ceil(frqpts/2));
            tv  = dt * ((1:pwpts) - ceil(pwpts/2));
            A   = zeros(frqpts,pwpts);
            
            for kfreq = 1:frqpts
                % generate a frequency offset pulse
                fop = p.frequencyOffsetPulse(fv(kfreq));
                [crossa,crossk] = p.computeCrossCorrelation(fop);
                
                [ipeak]=peakdet(abs(crossa),1.0e-6);
                crossa = abs(crossa(ipeak));
                crossk = crossk(ipeak);
                
                % convert crossk to time
                crossk = crossk * p.Ts;
                % pick off amplitude values from this data
                A(kfreq,:) = interp1(crossk,crossa,tv);
            end
            
            if ( nargout == 0 )
                %surf(1e6*tv,1e-3*fv,10*log10(abs(A)),'LineStyle','none');
                surf(1e6*tv,1e-3*fv,10*log10(A),'LineStyle','none');
                xlabel('Time Shift (\musec)');
                ylabel('Frequency Shift (kHz)');
                title('Correlation Magnitude (dB)');
                view(0,90);
                grid on;
                axis tight;
                colorbar;
            end
            
        end
        
        function SumPulse = plus(p,a)
            SumPulse = Pulse;
            SumPulse.amplitudeVector = p.amplitudeVector + a.amplitudeVector;
            SumPulse.timeVector = p.timeVector;
        end
        
    end
end