classdef FixedPoint
    
    properties (Access = public)
        value;
    end
    
    properties (Access = private)
        
        nIntBits;
        nFracBits;
        signBit;
        length;   % total number of bits
        
        % need to implement a way to specify overflow handling
        overflow; % could be "saturate" or "overflow";
        rounding = 'nearest'; % could be "floor", "ceiling", "nearest"
        
    end
    
    methods
        
        % goal is to have a class type defined and then can instantiate
        % based on this type
        % e.g.
        %   Q10_16 = FixedPoint(10,14,1);  % 25-bit number
        %   fpi = int16(pi)
        %   f = int
        %   y = 2*fpi*f
        
        function this = FixedPoint(nIntBits,nFracBits,signBit)
            % sign bit assumed to be included
            if(~exist('signBit'))
                signBit = 1;
            end
            this.nIntBits = nIntBits;
            this.nFracBits = nFracBits;
            
            this.length = nIntBits + nFracBits + signBit;
            
            this.value = [];
        end
        
        function obj = assign(this,val)
            
            % create new instance
            obj = this;
            % convert value to truncated fixed point format.
            intPart = floor(val);
            fracPart = mod(val,1);
            
            % truncate to fractional precision
            if(strcmp(obj.rounding,'floor'))
                
            elseif(strcmp(obj.rounding,'ceiling'))
                
            elseif(strcmp(obj.rounding,'nearest'))
                fixedPt = round(val * 2^obj.nFracBits) * 2^-obj.nFracBits;
            end
            
            % handle integer overflow
            if( fixedPt > (2^obj.nIntBits-1) ) % positive number
                if(strcmp(obj.overflow,'saturate'))
                    fixedPt = (2^obj.nIntBits-1) + mod(fixedPt,1);
                elseif(strcmp(obj.overflow,'overflow'))
                    % NEED TO DO...
                    if(obj.signBit == 1)
                        fixedPt = mod(fixedPt,2^obj.nIntBits-1)
                    else
                        
                    end
                end
            elseif( fixedPt < -2^obj.nIntBits ) % negative number
                if(strcmp(obj.overflow,'saturate'))
                    fixedPt = -2^obj.nIntBits - mod(fixedPt,1);
                elseif(strcmp(obj.overflow,'overflow'))
                    % NEED TO DO...
                    if(obj.signBit == 1)
                        fixedPt = mod(fixedPt,2^obj.nIntBits-1)
                    else
                        
                    end
                end
            end
            
            obj.value = fixedPt;
            
        end
        
        function obj = plus(fp1,fp2)
            
            % if the two fixed point numbers being added together are not
            % the same, then need to determine which format to output via
            % an optional argument (being the class object you want the
            % result to inherit from). If this optional argument is not
            % present, then the output result will be generated such that
            % the result will be accomodated.
            % BUT THIS IS TO BE AN OVERLOADED FUNCTION!!! HOW TO GET THIS
            % INFORMATION INTO THE GENERIC "fpsum = fp1 + fp2" FORMAT???
            
            if(class(outputFormat,'FixedPoint'))
                
            end
            
        end
        
        function val = get(this,str)
            val = this.(str);
        end
        
        function obj = convert(old,new)
            
        end
        
        
    end
    
    
    
end