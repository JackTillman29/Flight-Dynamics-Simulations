classdef LSFR < handle
    properties
        nBits;
        initVal;
        tapVal;
        polynomial;
        currentVal;
    end
    methods
        function obj = LSFR(nBits,polynomial,initVal)
            
            obj.tapVal = 0;
            polynomial = polynomial(end:-1:1);
            for k = 1 : length(polynomial)
                obj.tapVal = obj.tapVal + 2^(k-1) * polynomial(k);
            end
            
            obj.nBits = nBits;
            if(log2(obj.tapVal+1) > nBits)
                error('Feedback value exceeds number of bits!');
            end
  
            if(log2(initVal+1) > nBits)
                error('Initial value exceeds number of bits!');
            end
            obj.initVal = initVal;
            obj.currentVal = initVal;
        end
        function step(this)
            fbv = mod(bitand(this.currentVal,this.tapVal),2);
            fbv = bitshift(fbv,this.nBits-1,this.nBits);
            this.currentVal = bitshift(this.currentVal,-1) + fbv;
        end
        
        function display(this)
            fprintf(['Taps : ' dec2bin(this.tapVal,this.nBits) '\n']);
            fprintf(['Value: ' dec2bin(this.currentVal,this.nBits) '(' num2str(this.currentVal) ')\n']);
        end
    end
end