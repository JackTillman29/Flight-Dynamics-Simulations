classdef FirFilterC  < handle
%disp(['Use:']);
%disp([' applyFir(<fltobj>, <input sample>)']);
%disp([' applyFirData(<fltobj>, <input sample set>)']);
%disp([' applyFirDataFast(<fltobj>, <input sample set>)']);
%disp([' resetFir(<fltobj>)']);
%disp([' [input,output] = stepFir(<fltobj>,<#samples>,<step high sample>,<step low sample>)']);
%disp([' [freq,output] = bodeFir(<fltobj>,<N freq points>)']);
%disp([' plot(<fltobj>)']);    
    properties
        tapgains = [];
        tapgains0 = []; % initial tap gains
        tapsignals = 0;
        ndelays = 0;
        filename = '';
    end
    methods
        
        function this = FirFilterC(filename)
            if(exist('filename'))
                if(~strcmpi(filename,'null_hilbert'))
                    this.tapgains = load(filename);
                    this.tapgains0 = this.tapgains;
                    this.tapsignals = 0*this.tapgains;
                    this.ndelays = length(this.tapgains);
                    this.filename = filename;
                end
            else
                warning('Creating empty FIR filter!.');
            end
            
            %disp(['Use:']);
            %disp([' applyFir(<fltobj>, <input sample>)']);
            %disp([' applyFirData(<fltobj>, <input sample set>)']);
            %disp([' applyFirDataFast(<fltobj>, <input sample set>)']);
            %disp([' resetFir(<fltobj>)']);
            %disp([' [input,output] = stepFir(<fltobj>,<#samples>,<step high sample>,<step low sample>)']);
            %disp([' [freq,output] = bodeFir(<fltobj>,<N freq points>)']);
            %disp([' plot(<fltobj>)']);

        end
        
        % If you would like to use this class, but provide your own taps
        function this = OverrideTaps(this,taps,newname)
            if(size(taps,2) > size(taps,1))
                taps = taps';
            end
            this.tapgains = taps;
            this.ndelays = length(taps);
            this.tapsignals = zeros(this.ndelays,1);
            this.filename = newname;
        end
        
        % Added to enable rapid iteration of filters within a script using
        % the frequency domain window method. Note: for a low pass, you
        % must specify the bandwidth in quadrature (e.g. a 1MHz low pass
        % spans from -1MHz to +1MHz
        function this = AutogenFir(this,Fs,NFIR,bStart,bStop,window_fcn,name)
            
            % If no arguments are passed in, dump the syntax hint.
            if(nargin == 1)
                disp('Syntax: AutogenFir(this,Fs,NFIR,bStart,bStop,window_fcn,name)');
                return;
            end
            
            % If no name is passed in, use the default
            if(~exist('name','var'))
                name = 'Untitled FIR';
            end
            
            % If no window function is passed in (or is empty), use the default
            % rectangular window
            if(~exist('window_fcn'))
                window_fcn = @rectwin;
            elseif(isempty(window_fcn))
                window_fcn = @rectwin;
            end
            
            % To do this adequately for shorter FIR filters, need to
            % oversample in the frequency domain so as to achieve finer
            % frequency resolution in design phase.
            %  (NFFT >> NFIR)  or (NFFT,NFIR must be much greater than 0)
            oversamp_fac = 30;
            
            % Generate centered frequency vector from -Fs/2 to Fs/2
            df = Fs/(NFIR*oversamp_fac);
            fvec = [-(Fs/2):df:(Fs/2-df)];
            
            % Create the reference FFT window based upon input bandwidth
            desFft = double((fvec > bStart) & (fvec < bStop));
            
            if(sum(desFft) == 0)
                error('Too narrow of a filter requested!');
            end
            
            % shift into 0 to Fs
            desFft = fftshift(desFft);
            
            % IFFT the desired window and apply the user defined window
            % function to the filter impulse response (this will lower the
            % sidelobes of the filter in the frequency domain)
            desTap = ifftshift(ifft(desFft)).*window_fcn(NFIR*oversamp_fac);
            
            % Now window the taps to the requested length NFIR
            if(mod(NFIR,2) == 1)
                % if odd length NFIR
            %     desTap = ifftshift(desTap);
                idx_center = ceil(NFIR*oversamp_fac / 2);
                idx = [-(NFIR-1)/2:(NFIR-1)/2] + idx_center + 1;
            else
                % if even length NFIR
                idx_center = round(NFIR*oversamp_fac / 2);
                idx = [-(NFIR/2):((NFIR/2)-1)] + idx_center;
            end
            % USEFUL debugging
%                 ncen = [1:(NFIR*oversamp_fac)];
%                 figure;
%                 plot(ncen,desTap);
%                 hold on;
%                 plot(idx,desTap(idx),'r');
%                 
%                 figure; plot(fvec,fftshift(desFft))
%                 hold on;
%                 plot(fvec,abs(fftshift(fft(desTap,NFIR*oversamp_fac))))

            desTap = desTap(idx); % select only the windowed portion
            
            % Apppy the final impulse response to the filter object
            this.OverrideTaps(conj(desTap),name); 
        end
        
        function this = shiftFir(this,Fc,Fs)
            % quadrature shift, if real data is to be filtered, only use
            % the real part of the taps
            Ts = 1/Fs;
            for k = 1 : length(this.tapgains)
              this.tapgains(k) = this.tapgains(k) * exp(1i*2*pi*Fc*k*Ts);
            end
        end
        
        % jah - new 2017-09-21
        function this = unshiftFir(this)
            this.tapgains = this.tapgains0;
        end
        
        function this = normalizeGains(this)
            this.tapgains = this.tapgains ./ abs(sum(this.tapgains));
        end
        
        function y = applyFir(this,u)
            this.tapsignals(2:this.ndelays) = this.tapsignals(1:(this.ndelays-1));
            this.tapsignals(1) = u;
            y = sum(this.tapsignals .* this.tapgains);
        end
        
%         function y = applyFirData(this,u)
%             resetFir(this);
%             y = 0*u;
%             for k = 1:length(u)
%                 y(k) = applyFir(this,u(k));
%             end
%         end
        
        function y = applyFirData(this,u)
            resetFir(this);
            y = filter(this.tapgains,1,u);
        end
        
        function y = applyDecFirData(this,u,decFactor)
            % decFactor should be an integer factor of original filter
            % Implements a polyphase decimator for efficiency
            % note: output Fs will be Fs(in) / decFactor
            
            % build the 1D input vector to be reshaped for polyphase
            temp = [zeros(1,decFactor-1) u zeros(1,1)];
            
            % determine how much to truncate the 1D data so it will fit
            % nicely into a 2D reshaped array
            trunc_fac = floor((length(u)+decFactor)/decFactor)*decFactor;

            % create the 2D reshaped (decimated) polyphase input
            u_final = flipud(reshape( ...
                temp(1:trunc_fac), ...
                decFactor, ...
                floor((length(u)+decFactor)/decFactor)));
            
            % throw away last column
            u_final = u_final(:,1:end-1);
            
            % initialize the polyphase sum to zero
            y = 0;
            
            % apply individual polyphase filters in their respective data
            % sets
            for k = 1 : decFactor
                % sum the output
                y = y + filter(this.tapgains(k:decFactor:end),1,u_final(k,:));
            end
            
            % Done!!
        end
        
        
        function y = applyFirDataFast(this,u)
            warning('applyFirDataFast is deprecated. Use applyFirData method. It is much faster.');
            resetFir(this);
            flipFlag = 0;
            if(size(u,1) < size(u,2))
                u = u.';
                flipFlag = 1;
            end
            N = this.ndelays + length(u) - 1;
            fft_taps = fft(this.tapgains,N);
            fft_u    = fft(u,N);
            
            y = ifft(fft_taps .* fft_u);
            y = y(1:length(u));
            if(flipFlag)
                y = y.';
            end
        end
        
        function resetFir(this)
            this.tapsignals = 0.0 * this.tapsignals;
        end
        
        function [u,y] = stepFir(this,N,iUp,iDown)
            u = zeros(1,N);
            y = zeros(1,N);
            for k = 1:N
                if( (k >= iUp ) && (k < iDown) )
                    u(k) = 1.0;
                else
                    u(k) = 0.0;
                end
                y(k) = applyFir(this,u(k));
            end
            if(nargout == 0)
                figure;
                plot([u.' y.']);
                grid on;
                xlabel('Sample');
                ylabel('Magnitude');
                title(['step respone of ' this.filename],'Interpreter','none');
            end
        end
        
        function [f,y] = bodeFir(this,N,varargin)
            
            f = 0:(N-1);
            resetFir(this);
            y = fft(this.tapgains,N);
            fscale = 1;
            fstring = 'k\times Fs';
            if(nargin > 2)
                fscale = varargin{1};
                fstring = varargin{2};
            end
            
            if(nargout == 0)
                figure;
                subplot(2,1,1)
                plot(fscale*f./N,20*log10(abs(y)));
                grid on;
                ylabel('Magnitude 20 log_{10}');
                title(['Frequency Response of ' this.filename],'Interpreter','none');
                
                subplot(2,1,2)
                plot(fscale*f./N,angle(y));
                grid on;
                xlabel(['Frequency (' fstring ')']);
                ylabel('Phase [rad]');
                
                
%                 figure;
%                 term1 = fscale*((0:(3*N-1))-N)./N;
%                 term2 = 20*log10([abs(y)' abs(y)' abs(y)']);
%                 plot(term1,term2);
%                 grid on;
%                 xlabel(['Frequency (' fstring ')']);
%                 ylabel('Magnitude 20 log_{10}');
%                 title(['Frequency Response of ' this.filename],'Interpreter','none');
            end
        end
        
        function plot(this,N,varargin)
            try
                bodeFir(this,N,varargin{:});
            catch
                fprintf('Incorrect arguments: plot(this,N,varargin)\n');
            end
        end
        
        function varargout = plotTransient(this,Fs,inFreq,inPulse,inPri,nPulse)
            p = which('PWM_Signal');
            if(isempty(p))
                p = which('FirFilterC');
                try_path = strrep(p, ...
                    'DigitalFilters\Iowa Hills MATLAB\FirFilterC.m', ...
                    'Signal Generation\PWM Signal');
                try
                    addpath(try_path);
                catch
                    error('Unable to locate PWM_Signal! (Required for this method)');
                end
            end
            Ts = 1/Fs;
            N = floor(nPulse * inPri / Ts);
            
            u = PWM_Signal(inPulse,inPri,Ts,N,1);
            tIn = (0:length(u)-1) * Ts;
            modu = exp(1i*2*pi*inFreq*tIn);
            y = this.applyFirData(u.*modu);
            
            if(nargout == 1)
                varargout{1} = y;
            end
            
            figure;
            hp = plot( ...
                tIn,u,tIn,real(y),tIn,abs(y));
            grid on;
            xlabel('Time (sec)');
            ylabel('Amplitude');
            set(hp(1),'Color','b','LineWidth',2);
            set(hp(2),'Color',[1 1 1]*.7,'LineWidth',1);
            set(hp(3),'Color',[1 0 0]*.7,'LineWidth',1);
            legend('Input','Output(I)','Output(Env)','Location','best');
            title([this.filename ' - Time Domain Response']);
            
        end
        
        function plotZeros(this)
            z = roots(this.tapgains);
            figure
            hax = axes; hold(hax,'on'); axis equal;
            plot(real(z),imag(z),'ko');
            % plot unit circle
            th = linspace(0,2*pi,1000);
            plot(cos(th),sin(th),'k--');
        end
        
        function [mag,phase] = drive(this,kFreq)
            % evaluate steady state value for signal
            k = (0:(this.ndelays-1));
            
            sigout = sum(this.tapgains.' * exp(1i*2*pi*(kFreq'*k)).',1);
            mag = abs(sigout);
            phase = angle(sigout);
            
        end
        
        function writeOutputFile(this,filename)
           fid = fopen(filename,'w');
           fprintf(fid,'%18.6e\n',this.tapgains);
           fclose(fid);
        end
    end
end