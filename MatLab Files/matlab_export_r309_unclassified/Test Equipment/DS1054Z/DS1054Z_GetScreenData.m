function [out_data] = DS1054Z_GetScreenData(obj,fileprefix)
warning('off','instrument:fread:unsuccessfulRead');

try
    obj.InputBufferSize = 65535;
    obj.OutputBufferSize = 65535;
catch
end
if(~strcmp(obj.Status,'open'))
    warning('DS1054Z port NOT open! Attempting....');
    fopen(obj);
end

% these are oscope specific parameters
maxSamplesToRead = 250000;
headerSize = 12;


%% Get current state


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
fprintf(obj,':WAV:MODE NORM') % screen
fprintf(obj,':WAV:FORM ASC')

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
    
    
    fprintf(obj,':WAV:DATA?');
    [data,len]= fread(obj);
    
    real_data = sscanf(char(data(12:end)'),'%f,')';
    
    
    eval(['out_data.ch' num2str(ichan) ' = real_data;']);
    if(ichan == actChanNums(end))
    real_time = (0:length(real_data) - 1)*xinc + xorig;
    out_data.time = real_time;
    end
end

if(nargout == 0)
    disp('only works for ch1 right now');
    csvwrite([fileprefix '.csv'],[out_data.time' out_data.ch1']);
end

end

