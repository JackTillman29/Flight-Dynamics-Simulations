function [smpStart,smpStop,smpLen] = WaveDataPartition(totalSamples,blkSizeSmp)
% Determines how to break of sample data requests
sampleMemory = totalSamples; % set on device
if(blkSizeSmp >= totalSamples)
    blkSizeSmp = totalSamples;
end
maxReadAtOnce = blkSizeSmp;
nReads = ceil(sampleMemory/maxReadAtOnce);

startStopLenArray = zeros(nReads,3);

for kk = 1 : nReads
    if(kk == 1)
        startStopLenArray(1,:) = [1 maxReadAtOnce maxReadAtOnce];
        
    else
        startStopLenArray(kk,1) = startStopLenArray(kk-1,2) + 1;
        len = min([maxReadAtOnce sampleMemory - startStopLenArray(kk,1) + 1]);
        startStopLenArray(kk,2) = startStopLenArray(kk,1) + len - 1;
        startStopLenArray(kk,3) = len;
    end
end
smpStart = startStopLenArray(:,1);
smpStop  = startStopLenArray(:,2);
smpLen   = startStopLenArray(:,3);
end