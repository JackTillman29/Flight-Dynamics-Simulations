classdef WaveObj
    properties
        fft;
        frq;
        psd;
        cpsd;
        nfft;
        option;
        paths_good_flag = 0;
        name = '';
        tdStartTime;   % apply this offset to time arrays when plotting
        tdSignal;
        tdWindow;
        Fs;
        cumpow;
        iTimeStart;  % used for storing selected time start/stop on the STFT plot from interactive callback
        iTimeEnd;
        z0;  % characteristic impedance
        
        powerStr = 'dBW';
        powerStrFun = @(x) 10*log10(x);
        
        pwrUnitStr = 'dBW';
    end
    methods
        function this = WaveObj(tdSignal,tdWindow,Fs,NFFT,Option,varargin)
            % this = WaveObj(tdSignal,tdWindow,Fs,NFFT,Option,varargin)
            % defaults
            this.name = 'signal';
            this.tdStartTime = 0;
            this.z0 = 1;
            % variable inputs
            for i = 1:2:length(varargin)
                this.(varargin{i}) = varargin{i+1};
            end
            
            this = SetArray(this,tdSignal,'tdSignal');
            this = SetArray(this,tdWindow,'tdWindow');
            this.Fs = Fs;
            if(isempty(tdWindow))
                tdWindow = 0*tdSignal + 1;
            end
            if(~exist('NFFT','var'))
                NFFT = length(tdSignal);
            end
            if(~exist('Fs','var'))
                Fs = 1;
            end
            if(nargin == 5)
                [fft, frq, psd, cpsd] = WaveFft(tdSignal,tdWindow,Fs,NFFT,Option,'z0',this.z0);
                this.option = Option;
            else
                [fft, frq, psd, cpsd] = WaveFft(tdSignal,tdWindow,Fs,NFFT,[],'z0',this.z0);
                this.option = 'none';
            end
            this = SetArray(this,fft,'fft');
            this = SetArray(this,frq,'frq');
            this = SetArray(this,psd,'psd');
            this = SetArray(this,cpsd,'cpsd');
            this.nfft = NFFT;
            this.cumpow = this.cpsd(end);
        end
        
        function this = SetArray(this,in,field)
            % SetArray ensures that all arrays in this object are column
            % vectors.
            % check if array is row vector
            if(size(in,2) > 1)
                this.(field) = in.';
            else
                this.(field) = in;
            end
        end
            
        
        function this = SetName(this,name)
            this.name = name;
        end
        
        function this = SetZ0(this,z0)
            this.z0 = z0;
%             this.psd = this.psd / z0;
%             df = diff(this.frq,1);
%             df = df(1);
%             this.cpsd = cumsum(this.psd) .* df;
        end
        
        function this = SetStartTime(this,tdStartTime)
            this.tdStartTime = tdStartTime;
        end
        
        function varargout = plot(thisin,fSc,fU,fO,pwrmode,hax1,hax2)
            %  plot(thisin,fSc,fU,fO,pwrmode,hax1,hax2)
            %       fSc:  freq scale
            %       fU:   freq unit string
            %       fO:   freq offset [in units given by "fU"]
            %       hax1: input axes handle for magnitude plot
            %       hax2: input axes handle for phase plot
            
            show_power_in_legend = 0;
            if(~exist('fSc'))
                fSc = 1;
                fU  = 'Hz';
            end
            if(~exist('fO'))
                fO = 0;
            end
            if(~exist('pwrmode'))
                pwrmode = 'psd';
            end
                
            
%             for k = 1 : length(thisin)
%                 this = thisin(k);
%                 figure;
%                 subplot(2,1,1);
%                 plot( ...
%                     fSc*this.frq+fO,10*log10(this.psd));
%                 title([this.name 'Power Spectral Density']);
%                 xlabel(['Frequency (' fU ')']);
%                 ylabel('dBW/Hz');grid on;
%                 subplot(2,1,2);
%                 
%                 s1 = this.cpsd./this.cpsd(end);
%                 plot(fSc*this.frq+fO,s1);
%                 ha2 = gca;
%                 title([this.name 'Cumulative Power Spectra [' num2str(this.cpsd(end)) 'W]']);
%                 xlabel(['Frequency (' fU ')']);
%                 ylabel('% Distribution');
%                 grid on;
%                 add_menubar_print_callbacks;
%                 linkaxes( ...
%                     findobj('Type','axes','Parent',gcf,'-not','Tag','legend'), ...
%                     'x');
%                 
%                 % Define a context menu; it is not attached to anything
%                 hcmenu = uicontextmenu;
%                 
%                 % Define callbacks for context menu items that change linestyle
%                 % Define the context menu items and install their callbacks
%                 item1 = uimenu(hcmenu, 'Label', 'On-Screen Power', 'Callback', @get_on_screen_power);
%                 set(hcmenu,'UserData',{ha2,this.cumpow});
%                 set(ha2,'uicontextmenu',hcmenu);
%                 
%             end
            
%             if(length(thisin) > 1) % combined plot
                % compute colors
                legend_cell = {};
                legend_cell_psd = {};
                
                if(~exist('hax1'))
                    hfig = figure('Position',[520 308 560 420]);
                    %        this pos offset is to prevent addition of
                    %        toolbar from pushing top of window offscreen
                    hax1 = subplot(2,1,1);
                else
                    hfig = get(hax1,'Parent');
                end
                
                    
                %line_colors = hsv(length(thisin));
                for k = 1 : length(thisin)
                    if(strcmpi(pwrmode,'psd'))
                        hp=plot(hax1,fSc*thisin(k).frq+fO,10*log10(thisin(k).psd));
                    elseif(strcmpi(pwrmode,'pwr'))
                        hp=plot(hax1,fSc*thisin(k).frq+fO,10*log10(thisin(k).psd*thisin(k).Fs/thisin(k).nfft));
                    else
                        error('unknown mode');
                    end
                    if(k == 1)
                        line_colors = get(gca,'ColorOrder');
                    end
                    hold(hax1,'on')
%                     line_colors = get(hax1,'ColorOrder');
                    set(hp,'Color',line_colors(k,:));
                    if(show_power_in_legend == 1)
                        legend_cell = {legend_cell{:},[thisin(k).name ['  ' num2str(thisin(k).cpsd(end)) 'W']]};
                    else
                        legend_cell = {legend_cell{:},thisin(k).name};
                    end
                    legend_cell_psd = {legend_cell_psd{:}, thisin(k).name };
                end
                ha1c = hax1;
                hold(hax1,'off')
                if(strcmpi(pwrmode,'psd'))
                    title(hax1,'PSD');
                else
                    title(hax1,'Power');
                end
                grid(hax1,'on')
                if(length(thisin) > 1)
                    hlegend = legend(legend_cell_psd);
                    % define a callback to occur AFTER as zoom event has
                    % occurred.
                    hzoom = zoom(gcf);
                    set(hzoom,'ActionPostCallback',@(src,evt) cb_update_legend_with_power(src,evt,thisin,fSc))
                else
                    hzoom = zoom(gcf);
                    set(hzoom,'ActionPostCallback',@(src,evt) cb_update_title_with_power(src,evt,thisin,fSc))
                end
                
                if(~exist('hax2'))
                    hax2 = subplot(2,1,2);
                end
                for k = 1 : length(thisin)
                    % -  % distribution on y-axis
                    %s1 = 100 * thisin(k).cpsd./thisin(k).cpsd(end);
                    % dBm on y-axis
                    %s1 = 10*log10(thisin(k).cpsd/1000);
                    s1 = 100*thisin(k).cpsd./thisin(k).cpsd(end);
                    hp = plot(hax2,fSc*thisin(k).frq+fO,s1);
                    hold(hax2,'on');
                    set(hp,'Color',line_colors(k,:));
                end
                ha2c = hax2;
                hold(hax2,'off');
                if(length(thisin) == 1)
                    if(~isempty(thisin.name))
                        title([thisin.name ' - Cumulative Power Spectra [' num2str(10*log10(thisin.cpsd(end))+30) 'dBm]']);
                    else
                        title(['Cumulative Power Spectra [' num2str(10*log10(thisin.cpsd(end))+30) 'dBm (' ...
                            num2str(10*log10(thisin.cpsd(end))) 'dBW)']);
                    end
                    title([thisin.name ' - Cumulative Power Spectra [' num2str(thisin.cpsd(end)) 'W]']);
                else
                    title(hax2,'CPSD');
                end
                grid(hax2,'on')
                if(length(thisin) > 1)
                legend(legend_cell,'Location','NorthWest');
                end
                linkaxes([ha1c ha2c],'x');
%             end

            xlabel(hax1,['Frequency (' fU ')']);
            if(strcmpi(pwrmode,'psd'))
                ylabel(hax1,'dBW/Hz');
            else
                ylabel(hax1,'dBW');
            end
            
            xlabel(hax2,['Frequency (' fU ')']);
             ylabel(hax2,'% Distribution');
           % ylabel(hax2,'dBm')
            add_menubar_print_callbacks(hfig);
            add_analysis_callbacks(hfig);
            
            if(nargout == 2)
                varargout{1} = hax1;
                varargout{2} = hax2;
            end
            
            % save wave object(s) to figure user data
            userData = get(gcf,'UserData');
            % add a new field to the user data structure
            userData.waveObjs = thisin;
            % save to figure user data field
            set(gcf,'UserData',userData);
            
            function cb_update_legend_with_power(src,evt,thisarr,fSc)
                % get axes handle
                try
                    haxs = findobj('Parent',src,'Type','Axes');
                catch
                    % for backwards compatibility
                    haxs = findobj('Parent',gcf,'Type','Axes','-not','Tag','legend');
                end

                hax = haxs(1);
                
                % get freq start/stop indices
                freqAxis = 'x';
                [iFreqStart iFreqEnd] = getFreqStartEnd(thisarr(1),hax,fSc,freqAxis);
                
                legend_cell = {};
                for k = 1:length(thisarr)
                    cpsd_win(k) = thisarr(k).cpsd(iFreqEnd) - thisarr(k).cpsd(iFreqStart);
                    legend_cell = {legend_cell{:}, ...
                        [thisarr(k).name ['  ' num2str(thisarr(k).powerStrFun(cpsd_win(k))) ' ' thisarr(k).powerStr]]};
                end
                
                % create/update legend on each axes
                for k = 2:length(haxs)
                    legend(haxs(k),legend_cell);
                end
                
            end
            
            function cb_update_title_with_power(src,evt,thisarr,fSc)
                % get axes handle
                try
                    haxs = findobj('Parent',src,'Type','Axes');
                catch
                    % for backwards compatibility
                    haxs = findobj('Parent',gcf,'Type','Axes');
                end
                
                % find axes with CPSD in the title
                for k = 1:length(haxs)
                    titlestr = get(haxs(k),'Title');
                    if(ishandle(titlestr))
                        titlestr = get(titlestr,'String');
                    end
                    if(~isempty(strfind(titlestr,'Cumulative Power Spectra')))
                        hax = haxs(k);
                        break
                    end
                end
                
                % get freq start/stop indices
                freqAxis = 'x';
                [iFreqStart iFreqEnd] = getFreqStartEnd(thisarr(1),hax,fSc,freqAxis);
                
                % set title string to power in this frequency band
                cpsd_win(k) = thisarr(k).cpsd(iFreqEnd) - thisarr(k).cpsd(iFreqStart);
                if(~isempty(thisarr(k).name))
                    title(hax,[thisarr(k).name ' - Cumulative Power Spectra [' num2str(10*log10(cpsd_win(k))) 'dBW]']);
                else
                    title(hax,['Cumulative Power Spectra [' num2str(10*log10(cpsd_win(k))) 'dBW]']);
                end
            end
            
        end
        
        function power = getPower(this,fStart,fStop)
            % returns power within a frequency band based on integrated
            % power spectral density
            idx_lo = find(this.frq>fStart,1,'first');
            idx_hi  = find(this.frq<fStop,1,'last');
            power = this.cpsd(idx_hi)-this.cpsd(idx_lo);
        end
        
        function varargout = fast_stft(this,windowLen,windowFcn,NFFT,tSc,tU,fSc,fU,fO,shiftflag,hax)
            % fast_stft(this,windowLen,windowFcn,NFFT,tSc,tU,fSc,fU,fO,shiftflag,hax)
            
            if(nargin == 1)
                windowLen = 256;
                windowFcn = [];
                NFFT      = windowLen;
                tSc       = 1e6;
                tU        = '\mus';
                fSc       = 1e-6;
                fU        = 'MHz';
                fO        = 0;
                shiftflag = 'noshift';
            end
            
            % input switching
            if(~exist('tSc'))
                tSc = 1;
                tU  = 's';
            end
            if(~exist('fSc'))
                fSc = 1;
                fU  = 'Hz';
            end
            if(~exist('shiftflag'))
                shiftflag = 'noshift';
            end
            if(~exist('fO'))
                fO = 0;
            end
            
            if(isempty(NFFT))
                NFFT = windowLen;
            end
            
            if(NFFT < windowLen)
                warning('NFFT < window length! You are throwing away data!');
            end
            
            if(isempty(windowFcn))
                wnd = ones(1,windowLen);
                %error('empty window function is broken at the moment!!');
            else
                wnd = windowFcn(windowLen);
            end
            
            if(exist('hax'))
                axes(hax);
            end
            
            
            % First, get data repackage parameters (but don't create new
            % data)
            N = length(this.tdSignal);
            nRows = floor(N/windowLen);
            
            if(strcmpi(shiftflag,'shift'))
                sVec = (linspace(-0.5,0.5,windowLen) * this.Fs) * fSc;
            else
                sVec = (((0:windowLen-1)./(windowLen-1)) * this.Fs) * fSc - fO;
            end
            tVec = (windowLen * (0:nRows-1) ./ this.Fs) * tSc;
            
            nSamples = nRows * windowLen;
            
            indat = reshape(this.tdSignal(1:nSamples),windowLen,nRows);
            wnd = window2d(indat,wnd,1);
            
            %figure;
            if(strcmpi(shiftflag,'shift'))
                imData = 20*log10( ...
                    abs(fftshift(fft( ...
                    indat.*wnd, ...
                    [], ...
                    1)./windowLen,1)).');
            else
                imData = 20*log10( ...
                    abs(fft( ...
                    indat.*wnd, ...
                    [], ...
                    1)./windowLen).');
            end
            hi=imagesc(sVec,tVec,imData);
            if(fO == 0)
                xlabel(['Frequency (',fU,')']);
            else
                xlabel(sprintf('Frequency \\Delta rel %4.1f (%s)',-fO,fU));
            end
            ylabel(['Time (' tU ')']);
            if(exist('hax'))
                set(hi,'Parent',hax);
            end
            title(sprintf('Frequency \\Delta rel %4.1f (%s)',-fO,fU));
              
        end
        
        function varargout = stft(this,windowLen,windowFcn,overlap,NFFT,...
                tSc,tU,fSc,fU,fO,hax)
            % stft(...)
            % Input Arguments:
            %   1.) WaveObj
            %   2.) analysis window length  (in samples)
            %       how many input samples to use (per frame)
            %   3.) windowFcn @hamming, etc, or empty
            %   3.) analysis window overlap (fraction between 0 and 1)
            %   4.) NFFT (number of frequency bins)
            %       if NFFT > window length, zero padding will occur
            %   5.) tSc: time scale factor
            %   6.) tU:  time unit string
            %   7.) fSc: freq scale factor
            %   8.) fU:  freq unit string
            %   9.) fO:  freq offset (in Hz);
            %  10.) hax: axes handle to where you want to plot the STFT image.

            % DEFAULTS FOR WHEN YOU ONLY WANT TO INPUT THE WaveObj
            if(nargin == 1)
                default_num_frames = 20;
                windowLen = length(this.tdSignal)./default_num_frames;
                windowFcn = [];
                overlap   = 0.1;
                NFFT      = windowLen;
                tSc       = 1e6;
                tU        = '\mus';
                fSc       = 1e-6;
                fU        = 'MHz';
                fO        = 0;
            end
            
            % input switching
            if(~exist('tSc'))
                tSc = 1;
                tU  = 's';
            end
            if(~exist('fSc'))
                fSc = 1;
                fU  = 'Hz';
            end
            if(~exist('fO'))
                fO = 0;
            end
            
            if(isempty(NFFT))
                NFFT = windowLen;
            end
            
            if(NFFT < windowLen)
                warning('NFFT < window length! You are throwing away data!');
            end
            
            if(isempty(windowFcn))
                wnd = ones(1,windowLen);
            else
                wnd = windowFcn(windowLen);
            end
            
            x = this.tdSignal;
            w = this.tdWindow;
            if(size(x,1) > size(x,2))
                x = x.';
                w = w.';
            end
            Fs = this.Fs;
            dt = 1/Fs;
            N = length(x);
            overlapSamples = floor(overlap*windowLen);
%             numframes = ceil(N/(windowLen - overlapSamples))
            %               (1.0 - %)*W
            %      W     |----|
            %  [--------]
            %        [--------]
            %            ....
            % N = length of time history
            % % = percent overlap [val between 0-1]
            % W = window length of STFT frame (in samples)
            %
            % total length:
            %   W + N*(1.0 - %)*W = L
            numframes = ceil((N-windowLen)/windowLen/(1-overlap)) + 1;
            
            
            idx1 = 1;
            idx2 = windowLen;
            cm = bone(256);
            hw = waitbar(0,sprintf('Computing STFT Frame %d of %d',1,numframes));
            
            for i = 1:numframes
                cmi = floor(i/numframes * 255) + 1;
                if(idx2 <= N)
                    xsel = x(idx1:idx2);
                    
                else
                    idx2 = N;
%                     xsel = [x(idx1:idx2) zeros(1,1)];
                    xsel = x(idx1:idx2);
                    Nlast = windowLen - length(xsel);
                    xsel = [xsel zeros(1,Nlast)];
                end
                
               [xfft, frq, psd, cpsd] = ...
                   WaveFft(xsel,wnd,Fs,NFFT);
                
                savefft(i,1:NFFT) = xfft;
                idx1 = idx2 + 1 - overlapSamples;
                idx2 = idx1 + windowLen - 1;
                try
                waitbar(i/numframes,hw,sprintf('Computing STFT Frame %d of %d',i,numframes));
                end
%                 set(hw,'Color',cm(cmi,:));
            end
            try
            close(hw);
            end
            
            t = [0:(numframes-1)]./numframes;
            t = t * dt * length(x) + this.tdStartTime;
            sig = 20*log10(abs(savefft));
            if(~exist('hax'))
                figure('Position',[520 348 560 420]);

                hi=imagesc(fSc*frq+fO,tSc*t,sig);
            else
                % check for existing image
                ex_im=findobj('Type','Image','Parent',hax);
                if(isempty(ex_im))
                    hi=imagesc(fSc*frq+fO,tSc*t,sig);
                    set(hi,'Parent',hax);
                else
                    set(ex_im,'CData',sig);
                    return;
                end
            end
            
            if(nargout == 1)
                varargout{1} = gca;
            end
            % add colorbar axis adjustment callback
            set(gcf,'KeyPressFcn',@cb_modify_colorbar_caxis)
            
            % assign cb_stft_time_select callback to image
            % check to ensure a uicontextmenu doesn't already exist
            if(~isempty(get(hi,'UIContextMenu')))
                hcmenu = get(hi,'UIContextMenu');
            else
                hcmenu = uicontextmenu;
            end
            item1 = uimenu(hcmenu,'Label','FFT of waveform at this time',...
                'Callback',...
                @(src,evt)cb_stft_to_fft(src,evt,this,tSc,fSc,fU,fO));
            item2 = uimenu(hcmenu,'Label','Time-domain signal at this time',...
                'Callback',...
                @(src,evt)cb_fft_to_td(src,evt,this,tSc,tU,fSc,fU,fO));
            set(hi,'UIContextMenu',hcmenu);
            
%             colormap(gray)
%             colormap(bone)
%             colormap()
            sig(find(isinf(sig))) = -500;
            meanVal = mean(mean(sig));
            minVal = min(min(sig));
            maxVal = max(max(sig));
            caxis([meanVal maxVal]);
            hc=colorbar();
            %set(hc,'YColor','w');
            if(isreal(x))
                xlim([0+fO fSc*Fs/2+fO]);
            else
                xlim([0+fO fSc*Fs+fO]);
            end
            xlabel(['Frequency (' fU ')']);
            ylabel(['Time (' tU ')']);
            title('Short Time Fourier Transform');
            add_menubar_print_callbacks;
            add_analysis_callbacks;
            
            % also, add these callbacks to buttons that get put into the
            % toolbar
            htb = findobj('Parent',gcf,'Type','uitoolbar');
            if(~isempty(htb))
                % get directory where this function is stored (which should contain misc\)
                waveobj_dir = fileparts(mfilename('fullpath')); 
                item4 = uipushtool(htb,'CData',imread([waveobj_dir,'\misc\fft_to_time_domain_button.png']),'TooltipString','FFT --> Time-Domain',...
                    'ClickedCallback',@(src,evt)cb_fft_to_td(src,evt,this,tSc,tU,fSc,fU,fO));
                item5 = uipushtool(htb,'CData',imread([waveobj_dir,'\misc\stft_to_fft_button.png']),'TooltipString','STFT --> FFT',...
                    'ClickedCallback',@(src,evt) cb_stft_to_fft(src,evt,this,tSc,fSc,fU,fO));
            end
            
            function cb_stft_to_fft(src,evt,this,tSc,fSc,fU,fO)
                % get image handle
                him = findobj('Parent',gca,'Type','Image');
                
                % get axes handle
                hax = get(him,'Parent');
                
                % get time indices
                timeAxis = 'y';
                [iTimeStart iTimeEnd t] = getTimeStartEnd(this,hax,tSc,timeAxis);
                
                % plot
                if(size(this.tdSignal,1) == 1)
                    rectWin = ones(1,iTimeEnd-iTimeStart+1);
                else
                    rectWin = ones(iTimeEnd-iTimeStart+1,1);
                end
                tmp = WaveObj(this.tdSignal(iTimeStart:iTimeEnd),...
                              rectWin,this.Fs,length(rectWin));
                plot(tmp,fSc,fU,fO);
            end
            
            function cb_fft_to_td(src,evt,this,tSc,tU,fSc,fU,fO)
                % get image handle
                him = findobj('Parent',gca,'Type','Image');
                
                % get axes handle
                hax = get(him,'Parent');
                
                % get time indices
                timeAxis = 'y';
                [iTimeStart iTimeEnd t] = getTimeStartEnd(this,hax,tSc,timeAxis);
                
                % plot
                figure;
                haxnew = axes;
                tmpt = t(iTimeStart:iTimeEnd) * tSc;
                tmpsig = this.tdSignal(iTimeStart:iTimeEnd);
                if(isreal(this.tdSignal))
                    hp = plot(tmpt,tmpsig);
                else
                    hp = plot(tmpt,...
                        [real(tmpsig) imag(tmpsig)]);
                    set(hp,'Color',0.7*ones(1,3));
                    set(hp(2),'LineStyle','--','LineWidth',1);
                    hold on; grid on;
                    plot(tmpt,abs(tmpsig),'Color',[0.7 0 0],'LineWidth',2');
%                     plot(t,-abs(this.tdSignal),'k','LineWidth',2')
                end
                xlabel(['Time (',tU,')']);
                title(this.name);
                legend('I','Q','Envelope');
                % add a callback to this plot for plotting FFT of windowed
                % data in this figure
                if(~isempty(get(haxnew,'UIContextMenu')))
                    hcmenu = get(haxnew,'UIContextMenu');
                else
                    hcmenu = uicontextmenu;
                end
                item1 = uimenu(hcmenu,'Label','FFT of visible time data',...
                    'Callback',@(src,evt) cb_fft_of_visible_data(src,evt,this,tSc,tU,fSc,fU,fO));
                set(haxnew,'UIContextMenu',hcmenu);
                % also, add these callbacks to buttons that get put into the
                % toolbar
                add_menubar_print_callbacks(gcf);
                add_analysis_callbacks(gcf);
                add_analysis_callbacks(gcf);
%                 item2 = uipushtool(htb,'CData',imread('misc\accumulate_fft_button.png'),'TooltipString','Accumulate FFTs (not working)',...
%                     'ClickedCallback',@(src,evt) cb_fft_of_visible_data(src,evt,this,tSc,tU,fSc,fU,fO));
                
                
                
            end
            
            function cb_fft_of_visible_data(src,evt,this,tSc,tU,fSc,fU,fO)
                
                % get time indices
                timeAxis = 'y';
                [iTimeStart iTimeEnd t] = getTimeStartEnd(this,gca,tSc,timeAxis);
                
                % plot
                if(size(this.tdSignal,1) == 1)
                    rectWin = ones(1,iTimeEnd-iTimeStart+1);
                else
                    rectWin = ones(iTimeEnd-iTimeStart+1,1);
                end
                tmp = WaveObj(this.tdSignal(iTimeStart:iTimeEnd),...
                              rectWin,this.Fs,length(rectWin));
                plot(tmp,fSc,fU,fO);
                
            end
            
            
        end
        
        function td_to_fft_of_visible_data(hax)
        
        end
        
        function stft_to_fft_of_visible_data(hax)
            
        end
        

        
    end
    
    % define private methods for use internally by the class
    methods(Access = private)
        
        function [iStart iEnd t] = getTimeStartEnd(this,hax,tSc,timeAxis)
            
            if(~exist('timeAxis'))
                timeAxis = 'x';
            end
            
            % get time data + startTimeOffset
            L = length(this.tdSignal);
            t = [0:L-1]./this.Fs + this.tdStartTime;

            if(strcmpi(timeAxis,'x'))
                timeLim = get(hax,'XLim');
            else(strcmpi(timeAxis,'y'))
                timeLim = get(hax,'YLim');
            end
            % get time units and convert back to seconds
            timeLim = timeLim ./ tSc;

            % find this window in original data
            iStart = find(t >= timeLim(1),1,'first');
            iEnd   = find(t <= timeLim(2),1,'last');
        end
        
        function [iStart iEnd] = getFreqStartEnd(this,hax,fSc,freqAxis)
            
            if(~exist('timeAxis'))
                freqAxis = 'x';
            end
            
            % get freq data + freqOffset
            freq = this.frq;

            if(strcmpi(freqAxis,'x'))
                freqLim = get(hax,'XLim');
            else(strcmpi(freqAxis,'y'))
                freqLim = get(hax,'YLim');
            end
            % get freq units and convert back to Hz
            freqLim = freqLim ./ fSc;

            % find this window in original data
            iStart = find(freq >= freqLim(1),1,'first');
            iEnd   = find(freq <= freqLim(2),1,'last');
        end
        
    end
    
    
    
end
