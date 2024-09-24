classdef biquad
% Implements transposed direct form II biquad filter
%
% Y(z)   b0 z^2 + b1 z + b2 
%----- = ------------------- -> Y(k) = -a1*Y(k-1)-a2*Y(k-2) + b0*U(k) + b1*U(k-1) + b2*U(k-2)
% U(z)   1 z^2  + a1 z + a2
    properties
        B =  [];
        A =  [];
        A0 = [];
        B0 = [];
        
        delay = [0 0];
        
        poles = [];
        zeros = [];
    end
    methods
        
        function this = biquad(b0,b1,b2,a1,a2)
            this.B =  [b0 b1 b2];
            this.A =  [1  a1 a2];
            this.B0 =  [b0 b1 b2];
            this.A0 =  [1  a1 a2];
            this.delay = [0 0];
            this.poles = roots(this.A);
            this.zeros = roots(this.B);
        end
        
        function this = resetFilter(this)
            this.delay = [0 0];
        end
        
        function y = applyFilter(this,u)
            y = filter(this.B,this.A,u);
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
            
            zi = []; % initial condition of filter... init to empty means all taps are zero
            y = filter(this.B,this.A,u,zi,signalDim);
        end
        
%         function y = applyFilter(this,u)
%             if(length(u) == 1)
%                 warning('Input of size one!!! This routine loops over a vector of input samples!!');
%             end
%             
%             % CHECK TO SEE IF mex files ARE IN THE CURRENT PATH (AND
%             % COMPILED?)
%             % in args: inArray, Bgains, Agains
%             % out args: outArray
%             try
%                 if(isa(u,'double'))
%                     y = iirBiquadStageApplyFilterDouble(u,this.B,this.A);
%                 elseif(isa(u,'single'))
%                     y = iirBiquadStageApplyFilterSingle(u,this.B,this.A);
%                 end
%             catch
%                 warning('Mex files are non-functional. They should be in your path with biquad_stages. Using built-in code (slower)');
%                 n = length(u);
%                 y = 0*u;
%                 for k = 1 : n
%                     y(k) =          this.B(1) * u(k) + this.delay(1);
%                     this.delay(1) = this.B(2) * u(k) - this.A(2) * y(k) + this.delay(2);
%                     this.delay(2) = this.B(3) * u(k) - this.A(3) * y(k);
%                 end
%                
%             end
%             
%         end
        
        function this = shiftFilter(this,Fc,Fs)
            % clear out the filter delays
            this=this.resetFilter();
            
            Ts = 1/Fs;
            
            k = 1; sc2 = exp(1i*2*pi*Fc*k*Ts);
            k = 2; sc3 = exp(1i*2*pi*Fc*k*Ts);
            
            this.B(2) = this.B(2) * sc2;
            this.A(2) = this.A(2) * sc2;
            this.B(3) = this.B(3) * sc3;
            this.A(3) = this.A(3) * sc3;
            
            % recompute poles and zeros
            this.poles = roots(this.A);
            this.zeros = roots(this.B);
        end
        
%         function this = shiftFilter(this,fracFs)
%             sc1 = exp(1i*2*pi*fracFs);
%             sc2 = exp(1i*4*pi*fracFs);
% 
%             this.B(2) = this.B(2) * sc1;
%             this.A(2) = this.A(2) * sc1;
%             this.B(3) = this.B(3) * sc2;
%             this.A(3) = this.A(3) * sc2;
%         end
        
        function this = unshiftFilter(this)
            this.B = this.B0;
            this.A = this.A0;
            % recompute poles and zeros
            this.poles = roots(this.A);
            this.zeros = roots(this.B);
             % clear out the filter delays
            this=this.resetFilter();
        end
        
        function y = evalFilter(this,s,Ts)
            z = exp(1i.*s.*Ts);
            y = evalFilterZ(this,z);
        end
        
        function y = evalFilterZ(this,z)
            y = polyval(this.B,z) ./ polyval(this.A,z);
        end
        
        function plot(this,npts,Fs,FsUnit)
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
            end
            da = 2*pi*(0:npts)./npts-pi;
            z  = exp(1i*da);
            y = evalFilterZ(this,z);
            figure;
            plot(fscale*da./pi,20*log10(abs(y)));
            if(Fs == 0)
                xlabel('Fraction of Sample Rate');
            else
                xlabel(['Frequency (' FsUnit ')']);
            end
            ylabel('20log_{10} A');
            grid on;
            title('Bi-Quad Filter Response');
        end
    end
end