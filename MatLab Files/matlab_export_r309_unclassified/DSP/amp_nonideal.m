classdef amp_nonideal < handle
    
    properties
        P1dB_W = inf; % default to inf (turns off)
        P1dB_mW = inf;
        P1dB_dBW = inf;
        P1dB_dBm = inf;
        
        gain_linear = [];
        gain_dB = [];
        
        intermod_dist = [];
        
        noise_figure = [];
        bw = [];
        
    end
    
    
    methods
        
        function setGain(this,G,varargin)
            % this.setGain(6)            --> sets gain to 6 dB
            % this.setGain(100,'linear') --> sets gain to 20 dB
            % this.setGain(20,'dB')      --> sets gain to 20 dB
            % set the amplifier gain
            %   assumes value in dB
            %   
            if(length(varargin) > 0)
                switch lower(varargin{1})
                    case 'db'
                        this.gain_dB = G;
                        this.gain_linear = 10^(G/10);
                    case 'linear'
                        this.gain_linear = G;
                        this.gain_dB = 10*log10(G);
                end
            else
                this.gain_dB = G;
                this.gain_linear = 10.^(G/10);
            end
        end
        
        function setP1dB(this,P1dB,varargin)
            if(length(varargin) == 1)
                units = varargin{1};
            else
                % using default 1-dB compression point input units of dBm.
                units = 'dbm';
            end
            switch lower(units)
                case 'dbm'
                    this.P1dB_dBm = P1dB;
                    this.P1dB_dBW = P1dB - 30;
                    this.P1dB_mW  = 10^(this.P1dB_dBm/10);
                    this.P1dB_W   = 10^(this.P1dB_dBW/10);
                case 'dbw'
                    this.P1dB_dBW = P1dB;
                    this.P1dB_dBm = P1dB + 30;
                    this.P1dB_mW  = 10^(this.P1dB_dBm/10);
                    this.P1dB_W   = 10^(this.P1dB_dBW/10);
                case 'mw'
                    this.P1dB_mW  = P1dB;
                    this.P1dB_W   = P1dB * 1e-3;
                    this.P1dB_dBm = 10*log10(this.P1dB_mW);
                    this.P1dB_dBW = 10*log10(this.P1dB_W);
                case 'w'
                    this.P1dB_W   = P1dB;
                    this.P1dB_mW  = P1dB * 1e3;
                    this.P1dB_dBm = 10*log10(this.P1dB_mW);
                    this.P1dB_dBW = 10*log10(this.P1dB_W);
            end
        end
        
        function setNoiseFigure(this,NF,varargin)
            if(length(varargin) == 1)
                units = varargin{1};
            else
                % without units given, assume noise figure is input in dB
                units = 'db';
            end
            switch lower(units)
                case 'db'
                    this.noiseFigure_dB = P1dB;
                    this.P1dB_dBm = P1dB + 30;
                    this.P1dB_mW  = 10^(this.P1dB_dBm/10);
                    this.P1dB_W   = 10^(this.P1dB_dBW/10);
                case 'linear'
                    this.P1dB_W   = P1dB;
                    this.P1dB_mW  = P1dB * 1e3;
                    this.P1dB_dBm = 10*log10(this.P1dB_mW);
                    this.P1dB_dBW = 10*log10(this.P1dB_W);
            end
        end
        
        function y = applyGain(this,x,varargin)
            y = sqrt(this.gain_linear) .* x;
            % apply simple clipping
            y(y > sqrt(this.P1dB_W))  = sqrt(this.P1dB_W);
            y(y < -sqrt(this.P1dB_W)) = -sqrt(this.P1dB_W);
        end
    end
    
    
end