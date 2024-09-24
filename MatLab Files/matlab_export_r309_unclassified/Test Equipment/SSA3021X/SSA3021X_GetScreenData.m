function [real_data,real_frq]=SSA3021X_GetScreenData(SSA,fileprefix)
% SSA = visa('ni','USB0::0xF4EC::0xEE38::XXXXXXXXXX::INSTR');
warning('off','instrument:fread:unsuccessfulRead');

SSA.InputBufferSize = 65535;
SSA.OutputBufferSize = 65535;

if(~strcmp(SSA.Status,'open'))
    warning('SSA3021X port NOT open! Attempting....');
    fopen(SSA);
end

% get start freq
scpi_sequence(SSA,'SENS:FREQ:START?');
fStart = str2num(char(fread(SSA)'));
scpi_sequence(SSA,'SENS:FREQ:STOP?');
fStop = str2num(char(fread(SSA)'));

scpi_sequence(SSA,':TRAC:DATA? 1');
ydata=fread(SSA);
real_data = sscanf(char(ydata'),'%f,')';
real_frq  = linspace(fStart,fStop,length(real_data));

if(nargout == 0)
    csvwrite([fileprefix '.csv'],[real_frq' real_data']);
end

end

