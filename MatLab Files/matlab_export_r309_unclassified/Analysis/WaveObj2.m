classdef WaveObj2 < handle
    % A lightweight alternative to the original beast
    % WaveObj2(<time domain signal>,parm1,val1,parm2,val2,...)
            % Parameters    [default]
            %  z0           1 Ohm
            %  Fs           1 Hz
            %  tdWindow     @rectwin
    properties
        z0 = 1; % impedance cancels out 1 / (R) in power calculation such that
                                % P = V^2  (a useful simplification)
        Fs = 1;
        tdWindow = @rectwin;
        tdSignal = [];
        
        
    end
    
    methods
        
        function this = WaveObj2(tdSig,varargin)
            
            
            % Point to time domain waveform
            this.tdSignal = tdSig;
            
            if(size(tdSig,1) > size(tdSig,2))
                this.tdSignal = this.tdSignal.';
            end
            
            % Loop over options
            for k = 1 : 2 : length(varargin)
                set(this,varargin{k},varargin{k+1});             
            end
            
        end
        
        % Generic "setter" function
        function set(this,param,value)
            try
                eval(['this.' param ' = value;']);
            catch
                error(['Error setting parameter ' param '!']);
            end
        end
        
        function plot(this,varargin)
            % Syntax: plot(<waveobj2>,name1,val1,name2,val2,...)
            %
            % Options (*) indicates default : 
            %  'name'           'description for plots'
            %  'one sided'      1 = 0 to Fs/2 Data (e.g. real input)
            %                   0 = -Fs/2 to Fs/2 Data (complex or real)
            %                (*)    default based on data type (real = 1)
            %  'fft plot'       1 = produce an FFT plot(amplitude)
            %                (*)0 = no plot
            %  'power plot'  (*)1 = produce a spectrum analyzer plot(power)
            %                   0 = no plot
            %  'psd plot'    (*)1 = produce a PSD plot (power/Hz & CPSD)
            %                   0 = no plot
            %  'power dbw'      1 = Plot power in dBW
            %                (*)0 = Plot power in dBm
            %  'frequency scale'  = 'GHz' | 'MHz' | 'KHz' | (*)'Hz'
            %  'frequency reference' = DC Datum (for plotting)
            %                          Value is in 'frequency scale' units
            %                         (*) default is 0
            %  'voltage scale'    = (*)'V' | 'mV' | 'uV' | 'nV'
            %
            % Prototype:
            % plot(<this>,         ...
            %   'one sided',    0, ...
            %   'fft plot',     0, ...
            %   'power plot',   0, ...
            %   'psd plot',     0, ...
            %   'power dbw',    0, ...
            %   'frequency scale', 'Hz', ...
            %   'frequency reference', 0.0, ...
            %   'voltage scale', 'V');
            
            
            plotFlags.fft      = 0;
            
            % Set default sidedness to one for real signals (can be user
            % overridden)
            if(isreal(this.tdSignal))                
                plotFlags.oneSided = 1;
            else
                plotFlags.oneSided = 0;
            end
            plotFlags.power    = 1;
            plotFlags.psd      = 1;
            plotFlags.dBW      = 1;
            plotFlags.frqAxis  = 'hz';
            plotFlags.ampAxis  = 'v';
            name = '';
            fRef = 0;
            
            % User options
            for k = 1 : 2 : length(varargin)
                switch(lower(varargin{k}))
                    case 'one sided'
                        plotFlags.oneSided = varargin{k+1};
                    case 'fft plot'
                        plotFlags.fft = varargin{k+1};
                    case 'power plot'
                        plotFlags.power = varargin{k+1};
                    case 'psd plot'
                        plotFlags.psd = varargin{k+1};
                    case 'power dbw'
                        plotFlags.dBW = varargin{k+1};
                    case 'frequency scale'
                        plotFlags.frqAxis = lower(varargin{k+1});
                    case 'frequency reference'
                        fRef = varargin{k+1};
                    case 'voltage scale'
                        plotFlags.ampAxis = lower(varargin{k+1});
                    case 'name'
                        name = varargin{k+1};
                    otherwise
                        error(['Unknown option: ' varargin{k}]);
                end
            end
            
            % Map frequency checks
            switch(plotFlags.frqAxis)
                case 'ghz'
                    frqScale = 1e-9;
                    frqName  = 'Frequency (GHz)';
                case 'mhz'
                    frqScale = 1e-6;
                    frqName  = 'Frequency (MHz)';
                case 'khz'
                    frqScale = 1e-3;
                    frqName  = 'Frequency (KHz)';
                case 'hz'
                    frqScale = 1e-0;
                    frqName  = 'Frequency (Hz)';
                otherwise
                    error(['Unknown frequency option: ' plotFlags.frqAxis]);
            end
            
            
            % Map voltage checks
            switch(plotFlags.ampAxis)
                case 'v'
                    ampScale = 1e-0;
                    ampName  = 'Amplitude (V)';
                case 'mv'
                    ampScale = 1e3;
                    ampName  = 'Amplitude (mV)';
                case 'uv'
                    ampScale = 1e6;
                    ampName  = 'Amplitude (\muV)';
                case 'nv'
                    ampScale = 1e9;
                    ampName  = 'Amplitude (nV)';
                otherwise
                    error(['Unknown amplitude option ' plotFlags.ampAxis]);
            end
            
            % Check for window function
            try
                temp = this.tdWindow(1);
            catch
                addpath( ...
                    strrep( ...
                        which('WaveObj2'),'Analysis\WaveObj2.m','Windows'));
            end
            
            
            N = length(this.tdSignal);
            
            % Compute two-sided FFT & shift
            fftres = fftshift( ...
                fft(this.tdSignal .* this.tdWindow(N)));
            fftres = fftres ./ N;
            
            % Compute two-sided FFT Freq vector
            fftfrq = fftshift( ...
                this.Fs * (0:N-1) ./ N);
            fftfrq(fftfrq > (this.Fs/2)) = fftfrq(fftfrq > (this.Fs/2)) - this.Fs;
            
            % Apply user options
            if(plotFlags.oneSided == 1)
                idxKeep = find(fftfrq > 0);
                fftres = 2.0 * fftres(idxKeep); 
                fftfrq = fftfrq(idxKeep);
                pwrScale = 0.5 / (1 * this.z0);
            else
                pwrScale = 1 / (1 * this.z0);
            end
            
            if(~isreal(this.tdSignal))
                pwrScale = pwrScale / 2.0;
                if(plotFlags.oneSided)
                    warning('Complex input detected! One-Sided results are generally for real inputs only!');
                end
            end
            
            figs = [];
            
            if(plotFlags.fft)
                hf = figure; 
                ha = gca;
                plot(ha,frqScale*fftfrq-fRef,ampScale * abs(fftres));
                if(strcmp(name,''))
                    title('FFT');
                else
                    title(['FFT of ' name]);
                end
                xlabel(frqName);
                ylabel(ampName);
                grid on;
                figs = [figs hf];
            end
            
            
            if(this.z0 == 1)
                impString = '';
            else
                impString = sprintf('@ %d\\Omega',floor(this.z0));
            end
            
            if(plotFlags.dBW)
                pwrOffset = 0;
                pwrUnit   = sprintf('dBW %s',impString);
                psdUnit   = sprintf('dBW/Hz %s',impString);
            else
                pwrOffset = 30;
                pwrUnit   = sprintf('dBm %s',impString);
                psdUnit   = sprintf('dBm/Hz %s',impString);
            end
            
            if(plotFlags.power)
               hf = figure;
               ha = gca;
               
               hLine = plot(ha,frqScale*fftfrq+fRef,10*log10(pwrScale .* abs(fftres).^2)+pwrOffset);
               if(strcmp(name,''))
                    title('Spectrum Power');
                else
                    title(['Spectrum Power of ' name]);
                end
               xlabel(frqName);
               ylabel(pwrUnit);
               grid on;
               figs = [figs hf];
            end
            
            if(plotFlags.psd)
               hf = figure;             
               ha = subplot(2,1,1);
               
               hLine = plot(ha,frqScale*fftfrq+fRef,10*log10((N/this.Fs) * pwrScale .* abs(fftres).^2)+pwrOffset);
               if(strcmp(name,''))
                    title('PSD');
                else
                    title(['PSD of ' name]);
                end
               xlabel(frqName);
               ylabel(psdUnit);
               grid on; 
               
               hb = subplot(2,1,2);
               hLine2 = plot(hb,frqScale*fftfrq+fRef,cumsum(pwrScale .* abs(fftres).^2));
               if(strcmp(name,''))
                    title('CPSD');
                else
                    title(['CPSD of ' name]);
                end
               xlabel(frqName);
               ylabel(['W' pwrUnit(4:end)]);
               grid on; 
               figs = [figs hf];
               
               linkaxes([ha hb],'x');
               
               c = uicontextmenu;
               set(hLine,'uicontextmenu',c);
               m1 = uimenu(c,'Label','Occupied Channel Power','Callback', {@this.channelPowerFcn,hLine2,hLine});
               m2 = uimenu(c,'Label','Quick FIR Filter','Callback', {@this.quickFirFilter,hLine2,hLine});
               
            end
            
            for k = 1 : length(figs)
                figure(figs(k));
                try
                    add_print_callbacks;
                catch
                addpath( ...
                    strrep( ...
                        which('WaveObj2'),'Analysis\WaveObj2.m','Printing'));
                end
                add_analysis_callbacks;
            end
            
        end
        
        function channelPowerFcn(varargin)
            
            
            % Input Mapping
            hLine2 = varargin{4}; % CPSD
            hLine  = varargin{5}; % PSD
            hA_PSD  = get(hLine, 'Parent');
            hA_CPSD = get(hLine2,'Parent');
            xRef = get(hLine2,'XData');
            yRef = get(hLine2,'YData');
            
            gPts = ginput(2);
            
            x1 = gPts(1,1);
            x2 = gPts(2,1);
            
            % now find the nearest "actual" data points
            [~,idxClosest1] = min(abs(x1-xRef));
            [~,idxClosest2] = min(abs(x2-xRef));
            
            % user may have clicked right-to-left, so fix if needed
            if(idxClosest2 < idxClosest1)
                temp = idxClosest1;
                idxClosest1 = idxClosest2;
                idxClosest2 = temp;
            end
            
            % now subtract cumulative power at (2) from (1)
            channelPower = yRef(idxClosest2) - yRef(idxClosest1);
            channelPowerDBW = 10*log10(channelPower);
            channelPowerDBm = channelPowerDBW + 30;
            
            % Graphics Mapping
            psdLimY = ylim(hA_PSD);
            psdLimX = xlim(hA_PSD);
            cpsdLimY = ylim(hA_CPSD);
            cpsdLimX = xlim(hA_CPSD);
            
       
            
            hTB = annotation('textbox');
            set(hTB,'parent',hA_PSD,'Position', ...
                [xRef(idxClosest1) psdLimY(1) xRef(idxClosest2)-xRef(idxClosest1) diff(psdLimY)]);
            set(hTB,'Color',0.1*[1 1 1],'EdgeColor',0.5*[1 1 1],'FaceAlpha',0.25,'BackgroundColor',0.2*[1 1 1])
            
           % set(hTB,'Position',[idxClosest1
           str2 = sprintf('%4.1f dBW (%4.1f dBm)',channelPowerDBW,channelPowerDBm);
           set(hTB,'String',{'\bfOccupied','Channel','Power',str2},'VerticalAlignment','bottom');
           
           
 
        end
        
        function quickFirFilter(varargin)
            %================= duplicated  from above ===============
            % Input Mapping
            this = varargin{1};
            hLine2 = varargin{4}; % CPSD
            hLine  = varargin{5}; % PSD
            hA_PSD  = get(hLine, 'Parent');
            hA_CPSD = get(hLine2,'Parent');
            xRef = get(hLine2,'XData');
            yRef = get(hLine2,'YData');
            
            gPts = ginput(2);
            
            x1 = gPts(1,1);
            x2 = gPts(2,1);
            
            % now find the nearest "actual" data points
            [~,idxClosest1] = min(abs(x1-xRef));
            [~,idxClosest2] = min(abs(x2-xRef));
            
            % user may have clicked right-to-left, so fix if needed
            if(idxClosest2 < idxClosest1)
                temp = idxClosest1;
                idxClosest1 = idxClosest2;
                idxClosest2 = temp;
            end
            % ============= end duplicated ====================
            
            % Check for FirFilterC
            try
                qf = FirFilterC();
            catch
                addpath( ...
                    strrep( ...
                        which('WaveObj2'),'Analysis\WaveObj2.m','DigitalFilters\Iowa Hills MATLAB'));
            end
            qf = FirFilterC();
            qf.AutogenFir(this.Fs,length(xRef), ...
                idxClosest1 / length(xRef) * this.Fs - this.Fs/2, ...
                idxClosest2 / length(xRef) * this.Fs - this.Fs/2, ...
                @rectwin, ...
                'Quick FIR');
            this.tdSignal = qf.applyFirData(this.tdSignal);
            
            
        end
        
    end
    
    
    
    
end