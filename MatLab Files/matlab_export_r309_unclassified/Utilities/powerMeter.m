function datStruct = powerMeter(signal,signalName,datStruct,varargin)
% datStruct = powerMeter(signal,signalName,datStruct,varargin)
%    outputs a structure, "datStruct", that contains measured power
%    quantities of the input signal, "signal", stored in a field called
%    "signalName", which is itself a structure containing the numerical
%    values.


if(~exist('datStruct'))
    datStruct = [];
end
    
instPwrSig = abs(signal).^2;

if(~isfield(datStruct,signalName))
    datStruct.(signalName) = [];
    datStruct.(signalName).meanPwr_W = mean(instPwrSig);
    datStruct.(signalName).peakPwr_W = max(instPwrSig);
    datStruct.(signalName).meanPwr_dBW = 10*log10(datStruct.(signalName).meanPwr_W);
    datStruct.(signalName).peakPwr_dBW = 10*log10(datStruct.(signalName).peakPwr_W);
else
    % if the field already exists, then make this an array
    N = length(datStruct.(signalName));
    datStruct.(signalName)(N+1).meanPwr_W = mean(instPwrSig);
    datStruct.(signalName)(N+1).peakPwr_W = max(instPwrSig);
    datStruct.(signalName)(N+1).meanPwr_dBW = 10*log10(datStruct.(signalName)(N+1).meanPwr_W);
    datStruct.(signalName)(N+1).peakPwr_dBW = 10*log10(datStruct.(signalName)(N+1).peakPwr_W);
end




end