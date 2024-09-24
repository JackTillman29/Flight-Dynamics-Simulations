classdef biquad_stages
    % Implements transposed direct form II biquad filter
    %
    % Y(z)   b0 z^2 + b1 z + b2
    %----- = -------------------
    % U(z)   1 z^2  + a1 z + a2
    properties
        biquads = [];
        nstages = 0;
        ihfilename = [];
    end
    methods
        
        function this = biquad_stages(init_stages)
            if(isnumeric(init_stages))
                if(init_stages == 0)
                    error('This method is no longer supported. Please pass in a filename or # or stages');
                end
                this.nstages = init_stages;
            elseif(isstr(init_stages))
                this.ihfilename = init_stages;
                if( exist(init_stages) ~= 2 )
                    error(['Can not locate file: ' init_stages]);
                end
                this = this.read_file();
            else
                error('biquad_stages must be initialized with either 1) # of stages or 2) iowa hills filename');
            end
            
            %             addpath('mex');
            
        end
        
%         function this = iirFilterC(filenam)
%             this = biquad_stages(0);
%             this.ihfilename = filenam;
%             [a0,a1,a2,b0,b1,b2]=read_iowa_hills(filenam);
%             this.biquads = [];
%             for k = 1 : length(a0)
%                 this = this.add_stage( biquad(b0(k),b1(k),b2(k),a1(k),a2(k)));
%             end
%             this.nstages = length(a0);
%         end
        
        function this = read_file(this)
            [a0,a1,a2,b0,b1,b2]=this.read_iowa_hills_internal();
            this.biquads = [];
            for k = 1 : length(a0)
                this = this.add_stage( biquad(b0(k),b1(k),b2(k),a1(k),a2(k)));
            end
            this.nstages = length(a0);
        end
        
        function this = add_stage(this,bstage)
            this.biquads = [this.biquads bstage];
        end
        
        function this = shiftStages(this,Fc,Fs)
            % this = shiftStages(this,Fc,Fs)
            if(nargin ~= 3)
                error('check # of inputs:  shiftStages(this,Fc,Fs)');
                return;
            end
            for k = 1:this.nstages
                this.biquads(k) = shiftFilter(this.biquads(k),Fc,Fs);
            end
        end
        %         function this = shiftStages(this,fracFs)
        %             for k = 1:this.nstages
        %                 this.biquads(k) = shiftFilter(this.biquads(k),fracFs);
        %             end
        %         end
        
        function this = resetFilter(this)
            for k = 1:length(this.biquads)
                this.biquads(k).delay = [0 0];
            end
        end
        
        function this = unshiftStages(this)
            for k = 1:this.nstages
                this.biquads(k) = unshiftFilter(this.biquads(k));
            end
        end
        
%         function y = applyFilter(this,u)
%             n = length(u);
%             y = u;
%             for k = this.nstages : -1 : 1
%                 y = applyFilter(this.biquads(k),y);
%             end
%         end
        
        function y = applyFilter(this,u)
            y = u;
            for k = this.nstages : -1 : 1
                y = applyFilter(this.biquads(k),y);
            end
        end
        
        function y = applyFilterMultipleSignals(this,u,signalDim)
            % y = applyFilterMultipleSignals(this,u,signalDim)
            % HELP:
            % u can be a matrix, where "signalDim" should be the dimension
            %     along which the signal time histories lie.
            % Example:
            %   u is a 4x10 matrix, containing 4 signals, each of length-10
            %   therefore, signalDim should be set to 2 for the filter to
            %   be applied to all 4 signals simultaneously by MATLAB's
            %   filter command:
            %        y = filtObj.applyFilterMultipleSignals(u,2);
            y = u;
            for k = this.nstages : -1 : 1
                y = applyFilterMultipleSignals(this.biquads(k),y,signalDim);
            end
        end
        
        function y = evalFilter(this,s,Ts)
            y = 1;
            for k = 1 : this.nstages
                y = y .* this.biquads(k).evalFilter(s,Ts);
            end
        end
        
        function y = evalFilterZ(this,z)
            y = 1;
            for k = 1 : this.nstages
                y = y .* this.biquads(k).evalFilterZ(z);
            end
        end
             
        function plot(this,npts,Fs,FsUnit,band)
            % plot(this,npts,Fs,FsUnit,band)
            %     all inputs are OPTIONAL
            %     may need npts >> 1e4, but avoid 1e12...
            % band allows for npts to be spread over just frequencies in
            % the specified band [f1 f2] in Hz
            if(nargin < 2)
                disp(['Helper: plot(this,npts,Fs,FsUnit,band)']);
                disp(['Helper: plot(this,1000,1e6,''MHz'',[-1 1])']);
            end
            if(~exist('npts','var'))
                npts = 10001;
            else
                if(mod(npts,2) == 0)
                    npts = npts + 1;
                end
            end
            if(~exist('FsUnit','var'))
                FsUnit = 'Hz';
            end
            if(~exist('Fs','var'))
                Fs = 0;
                fscale = 0.5;
            else
                fscale = 0.5*Fs;
            end
            
            if(~exist('band','var'))
                %da = 2*pi*(0:npts)./npts-pi;
                da = linspace(-pi,pi,npts);
            else
                da = linspace(pi * band(1)/(Fs/2),pi * band(2)/(Fs/2),npts);
            end
            z  = exp(1i*da);
            figure;
            ah1=subplot(2,1,1);
            if(Fs == 0)
                xlabel('Fraction of Sample Rate');
            else
                xlabel(['Frequency (' FsUnit ')']);
            end
            
            ylabel('20log_{10} A');
            grid on;
            title('Bi-Quad Filter Amplitude Response');
            ah2=subplot(2,1,2);
            
            if(Fs == 0)
                xlabel('Fraction of Sample Rate');
            else
                xlabel(['Frequency (' FsUnit ')']);
            end
            ylabel('Phase (deg)');
            grid on;
            title('Bi-Quad Filter Phase Response');
            
            co = get(gca,'ColorOrder');
            
            for kf = 1 : length(this)
            y = evalFilterZ(this(kf),z);
            axes(ah1);hold on;
            plot(fscale*da./pi,20*log10(abs(y)),'Color',co(kf,:));
            hold off;
            axes(ah2);hold on;
            plot(fscale*da./pi,180/pi * angle(y),'Color',co(kf,:));
            end
            linkaxes([ah1 ah2],'x');
        end
        
        function step(this,N,Ts,TsUnit)
            disp('helper: step(<obj>,N,Ts,TsUnit)');
            if(~exist('Ts','var'))
                Ts = 1;
            end
            if(~exist('TsUnit','var'))
                TsUnit = 'Sample #';
            end
            u = ones(1,N);
            u(1) = 0;
            y = this.applyFilter(u);
            figure;
            plot(Ts*(0:N-1),[u.' y.']);
            xlabel(['Time (' TsUnit ')']);
            ylabel('Output');
            title([this.ihfilename ' step response'],'Interpreter','none');
            
            
            
        end
        
        function plotband(this,f1,f2,npts,Fs,FsUnit,opts)
            warning('may deprecate?? consider migrating to ''band'' option of plot...');
            % plot(this,f1,f2,npts,Fs,FsUnit,opts)
            %     f1, f2 in Hz IF Fs is supplied (-Fs/2,Fs/2) 
            %     otherwise, normalized on interval (-1:1)
            %     all inputs except f1 and f2 are OPTIONAL
            %     may need npts >> 1e4, but avoid 1e12...
            if(~exist('npts','var'))
                npts = 10000;
            end
            if(~exist('FsUnit','var'))
                FsUnit = 'Hz';
            end
            if(~exist('Fs','var'))
                Fs = 0;
                fscale = 0.5;
            else
                fscale = 0.5*Fs;
                f = linspace(f1,f2,npts);
            end
            if(~exist('opts'))
                opts = 'mag+phs';
            end
            % omega = 2 * pi * f = 2 * pi * f/Fs
            %     (fscale = Fs/2... 2 / Fs moved to omega, so divide)
            w1 = pi * f1 / fscale;
            w2 = pi * f2 / fscale;
            da = linspace(w1,w2,npts);
%             da = 2*pi*(0:npts)./npts-pi;
            z  = exp(1i*da);
            y = evalFilterZ(this,z);
            figure;
            if(strcmpi(opts,'mag+phs'))
                ah1=subplot(2,1,1);
            end
            if(strcmpi(opts,'mag+phs') || strcmpi(opts,'mag'))
                plot(fscale*da./pi,20*log10(abs(y)));
                if(Fs == 0)
                    xlabel('Fraction of Sample Rate');
                else
                    xlabel(['Frequency (' FsUnit ')']);
                end
                ylabel('20log_{10} A');
                grid on;
                title('Bi-Quad Filter Amplitude Response');
            end
            
            if(strcmpi(opts,'mag+phs'))
                ah2=subplot(2,1,2);
            end
            if(strcmpi(opts,'mag+phs') || strcmpi(opts,'phs'))
                plot(fscale*da./pi,180/pi * angle(y));
                if(Fs == 0)
                    xlabel('Fraction of Sample Rate');
                else
                    xlabel(['Frequency (' FsUnit ')']);
                end
                ylabel('Phase (deg)');
                grid on;
                title('Bi-Quad Filter Phase Response');
            end
            
            if(strcmpi(opts,'mag+phs'))
                linkaxes([ah1 ah2],'x');
            end
            
        end
        
        function plotPoleZero(this)
            figure;
            hax = axes; hold(hax,'on'); axis equal;
            for k = 1:this.nstages
                plot(real(this.biquads(k).zeros(1:2)),imag(this.biquads(k).zeros(1:2)),'ko')
                plot(real(this.biquads(k).poles(1:2)),imag(this.biquads(k).poles(1:2)),'kx')
            end
            % plot unit circle
            th = linspace(0,2*pi,1000);
            plot(cos(th),sin(th),'k--')
        end
        
        function bodeIir(this,N)
            % KDS-what is this providing vs plot?
            % evaluate freq response via poles and zeros
            resetFilter(this)
            magZ = zeros(this.nstages,N);
            magP = zeros(this.nstages,N);
            phsZ = zeros(this.nstages,N);
            phsP = zeros(this.nstages,N);
            %             n = linspace(-1,1,N); % full spectrum
            n = linspace(0,1,N); % half spectrum
            exptest = exp(1i*pi*n);
            for k = 1:this.nstages
                z1 = exptest-this.biquads(k).zeros(1);
                z2 = exptest-this.biquads(k).zeros(2);
                p1 = exptest-this.biquads(k).poles(1);
                p2 = exptest-this.biquads(k).poles(2);
                magZ(k,:) = abs(z1) .* abs(z2);
                magP(k,:) = abs(p1) .* abs(p2);
                phsZ(k,:) = angle(z1) + angle(z2);
                phsP(k,:) = angle(p1) + angle(p2);
            end
            bodeMag = prod(magZ,1) ./ prod(magP,1);
            bodePhs = sum(phsZ,1) -  sum(phsP,1);
            normFactor = max(bodeMag);
            bodeMag = bodeMag ./ normFactor;
            figure;
            subplot(2,1,1)
            plot(n,bodeMag)
            ylabel('|H(z)|')
            subplot(2,1,2)
            plot(n,bodePhs)
            ylabel('\angle H(z)')
            xlabel('Normalized Frequency')
        end
        
        function [a0,a1,a2,b0,b1,b2] = read_iowa_hills_internal(this)
            fid = fopen(this.ihfilename);
            s  = [];
            a0 = [];
            a1 = [];
            a2 = [];
            b0 = [];
            b1 = [];
            b2 = [];
            fgetl(fid);
            fgetl(fid);
            s=fgetl(fid);
            while(~strcmp(s,'Nth Order'))
                %disp(s);
                
                c1 = fgetl(fid);
                a0 = [a0 str2num(c1(3:end))];
                c1 = fgetl(fid);
                a1 = [a1 str2num(c1(3:end))];
                c1 = fgetl(fid);
                a2 = [a2 str2num(c1(3:end))];
                
                c1 = fgetl(fid);
                b0 = [b0 str2num(c1(3:end))];
                c1 = fgetl(fid);
                b1 = [b1 str2num(c1(3:end))];
                c1 = fgetl(fid);
                b2 = [b2 str2num(c1(3:end))];
                
                s=fgetl(fid);
                s=fgetl(fid);
                if(strcmp(s,''))
                    break;
                end
            end
            fclose(fid);
        end
        
        
        % MATLAB-only functionality
        %   can cascade multiple biquad_stages together into a single
        %   object. Nice for combining multiple filters into a single
        %   filter.
        function this = combine_biquads(this,old)
            % combine old biquad_stages filter with this biquad_stages
            % filter
            
            for k = 1:old.nstages
                this = this.add_stage(old.biquads(k));
            end
            this.nstages = length(this.biquads);
            
        
        end
        
        % Added this as an overloaded operator
        function newobj = plus(a,b)
            newobj = combine_biquads(a,b);
            newobj.ihfilename = 'combined';
        end
        
    end
end