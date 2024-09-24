function outputStruct = PulseCompressionXCORRsStruct(input,kernel,timealign,varargin)
% function [inputCrossCorr,kleadlag,kernelAutoCorr,fftInput,fftKernel,fftFreqVect] = PulseCompressionXCORRs(input,kernel,timealign='align' or 'no_align',varargin)
% varargin{1} = output unit scale factor of x axis. If omitted, it is considered
% "samples"

if(~exist('timealign'))
    timealign = 'align';
end

[inputCrossCorr,inputCrossCorrTrim,kleadlag,kernelAutoCorr,fftInput,fftKernel,fftFreqVect] = PulseCompressionXCORRs(input,kernel,timealign,varargin);

outputStruct.signalCrossCor = inputCrossCorr;
outputStruct.signalCrossCorTrim = inputCrossCorrTrim;
outputStruct.kleadlag = kleadlag;
outputStruct.kernalAutoCor = kernelAutoCorr;
outputStruct.signalFft = fftInput;
outputStruct.kernalFft = fftKernel;
outputStruct.frq = fftFreqVect;


end