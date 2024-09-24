classdef runningAvg
    properties
        counter = 0;
        Ex = [];
        Ex2 = [];
        Varx = [];
        Minx = []
        Maxx = [];
        setup = 1;
    end
    methods
        function this = runningAvg(varargin)
%             if(~exist('setup'))
%                 % setup(1): compute mean
%                 % setup(2): compute var
%                 this.setup = [1];
%             else
%                 this.setup = setup;
%             end
        end
        
        function this = newdata(this,x)
            
            this.counter = this.counter + 1;
            if(this.counter == 1)
                this.Ex = x;
                this.Ex2 = x.^2;
                this.Varx = this.Ex2 - this.Ex.^2;
                this.Minx = x;
                this.Maxx = x;
            else
                % average
                this.Ex = (this.counter-1)./this.counter * this.Ex ...
                    + 1./this.counter * x;
                % variance
%                 if(this.setup(1) == 1)
                this.Ex2 =(this.counter-1)./this.counter * this.Ex2 ...
                    + 1./this.counter * x.^2;
                this.Varx = this.Ex2 - this.Ex.^2;
                
                this.Minx = min(this.Minx, x);
                this.Maxx = max(this.Maxx, x);
                
%                 end
            end
            
        end
        
        function x = getdata(this)
            x = this.Ex;
        end
        
        function x = getex(this)
            x = this.Ex;
        end
        
        function x = getvarx(this)
            x = this.Varx;
        end
        
        function x = getstdx(this)
            x = sqrt(this.Varx);
        end
        
        function x = getminx(this)
            x = this.Minx;
        end
        
        function x = getmaxx(this)
            x = this.Maxx;
        end
        
        
        
        function this = reset(this)
            this.Ex = [];
            this.Ex2 = [];
            this.Varx = [];
            this.counter = 0;
            disp('running average object is reset.')
        end
        
    end
    
end