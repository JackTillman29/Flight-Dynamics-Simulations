function data_object = DS1054Z_FullCapture(obj,setDepthAmount)
%% This function requires that the SCPI I/F be "closed"
fclose(obj);
warning('off','instrument:fread:unsuccessfulRead');
depthOptions = 24e6 * [1 1/2 1/20 1/200 1/2000];
flag.setDepth = ~isempty(setDepthAmount);

if(flag.setDepth)
    if(sum(setDepthAmount == depthOptions))
        % good to go
    else
        depthOptions'
        error(['Invalid memory depth provided: ' num2str(setDepthAmount)]);
    end
end


% these are oscope specific parameters
maxSamplesToRead = 250000;
headerSize = 12;


%% Get current state
fopen(obj)

% ======================= channel configuration
fprintf(obj,':CHAN1:DISP?')
ch1.active = str2num(char(fread(obj)));
fprintf(obj,':CHAN2:DISP?')
ch2.active = str2num(char(fread(obj)));
fprintf(obj,':CHAN3:DISP?')
ch3.active = str2num(char(fread(obj)));
fprintf(obj,':CHAN4:DISP?')
ch4.active = str2num(char(fread(obj)));

channelNum = [1 2 3 4];
channelActive = [ch1.active  ch2.active  ch3.active  ch4.active];
nChannels = sum(channelActive);
actChanNums = find(channelActive ==1);
smpDivisorArray = [1 2 4 4];
smpDivisor = smpDivisorArray(nChannels);

% === memory

if(flag.setDepth == 0)
    fprintf(obj,':ACQ:MDEP?')
    mdepth = fread(obj);
    if(strcmpi(char(mdepth(1:4)'),'auto'))
        fclose(obj);
        error('AUTO mode is not supported');
    end
    
    mdepth = str2num(char(mdepth)');
else
    fprintf(obj,[':ACQ:MDEP ' num2str(setDepthAmount/smpDivisor)])
    mdepth = setDepthAmount/smpDivisor;
    
end





fclose(obj);

% Modify buffer size
bufferSize = min([maxSamplesToRead mdepth]) + headerSize;
obj.InputBufferSize = bufferSize;
fopen(obj);

%%


sampleMemory = mdepth; %1.2e6;
[smpStart,smpStop,smpLen] = WaveDataPartition(sampleMemory,maxSamplesToRead);

if(sampleMemory > 1e6)
    sMstring = [num2str(sampleMemory/1e6) 'M'];
elseif(sampleMemory > 1e3)
    sMstring = [num2str(sampleMemory/1e3) 'K'];
else
    sMstring = num2str(sampleMemory);
end

hw = waitbar(0,['Reading ' sMstring ' samples from DS1054Z CH' num2str(actChanNums(1))]);

%
fprintf(obj,':STOP')
fprintf(obj,':SYST:BEEP 0');

for ichan = actChanNums
    
    fprintf(obj,[':WAV:SOUR CHAN' num2str(ichan)])
    
    fprintf(obj,':WAV:XINC?');
    xinc = str2num(char(fread(obj)'));
    fprintf(obj,':WAV:XREF?');
    xref = str2num(char(fread(obj)'));
    fprintf(obj,':WAV:XOR?');
    xorig = str2num(char(fread(obj)'));
    
    fprintf(obj,':WAV:YINC?');
    yinc = str2num(char(fread(obj)'));
    fprintf(obj,':WAV:YREF?');
    yref = str2num(char(fread(obj)'));
    fprintf(obj,':WAV:YOR?');
    yorig = str2num(char(fread(obj)'));
    
    fprintf(obj,':WAV:MODE RAW')
    fprintf(obj,':WAV:FORM BYTE')
    fullData = zeros(1,sampleMemory,'single');
    %tic
    for iRequest = 1 : length(smpStart)
        startSample = smpStart(iRequest);
        stopSample = smpStop(iRequest);
        %fprintf('Requesting sample %d through %d\n',startSample,stopSample);
        fprintf(obj,[':WAV:STAR ' num2str(startSample)])
        fprintf(obj,[':WAV:STOP ' num2str(stopSample)])
        fprintf(obj,':WAV:DATA?');
        [data,len]= fread(obj,bufferSize);
        
        
        %fullData(startSample:stopSample) = (single(data(12:end-1)) - yref)*yinc;
        fullData(startSample:stopSample) = data(12:end-1);
        %waitbar(,hw);
        waitbar( ...
            iRequest./length(smpStart), ...
            hw, ...
            ['Reading ' sMstring ' samples from DS1054Z CH' num2str(ichan)]);
    end
    eval(['data_object.ch' num2str(ichan) '.xReference = xref;']);
    eval(['data_object.ch' num2str(ichan) '.yReference = yref;']);
    eval(['data_object.ch' num2str(ichan) '.xIncrement = xinc;']);
    eval(['data_object.ch' num2str(ichan) '.yIncrement = yinc;']);
    eval(['data_object.ch' num2str(ichan) '.xOrigin = xorig;']);
    eval(['data_object.ch' num2str(ichan) '.yOrigin = yorig;']);
    eval(['data_object.ch' num2str(ichan) '.rawData = fullData;']);
end
close(hw);

%toc
fprintf(obj,':SYST:BEEP 1');
fprintf(obj,':RUN');
fclose(obj)
end