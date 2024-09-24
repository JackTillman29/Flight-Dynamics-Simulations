classdef HilbertFilter
    properties
        nShiftReg
        iChanFIR % simple delay line
        qChanFIR % Hilbert FIR filter
    end
    methods
        function this = HilbertFilter(filename)
             this.qChanFIR = FirFilterC(filename);
             if(~mod(this.qChanFIR.ndelays,2))
                 error('Hilbert must be odd # of taps!');
             end
             
             % Iowa Hills code seems to generate a phase lead (vs lag)
             this.qChanFIR.tapgains = -1.0 * this.qChanFIR.tapgains;
             
             this.nShiftReg = (this.qChanFIR.ndelays+1)/2;
             
             this.iChanFIR = FirFilterC('null_hilbert');
             
             temp = zeros(1,this.nShiftReg);
             temp(end) = 1.0;
             
             OverrideTaps(this.iChanFIR,temp,'ichannel');
             
             this.iChanFIR.ndelays = this.nShiftReg;
             this.iChanFIR.tapsignals = zeros(this.nShiftReg,1);  
        end
        
        function [yi,yq] = applyHilbertData(this,u)
            yi = this.iChanFIR.applyFirData(u);
            yq = this.qChanFIR.applyFirData(u);
        end
        
        function [yi,yq] = applyHilbertDataFast(this,u)
            yi = this.iChanFIR.applyFirDataFast(u);
            yq = this.qChanFIR.applyFirDataFast(u);
        end
        
    end
end