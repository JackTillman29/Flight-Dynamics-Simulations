function SDG805_OnOff(SDG805,mode)
%SDG805 = visa( ...
%    'ni', ...
%    'USB0::0xF4ED::0xEE3A::SDG00003140697::INSTR');
%fclose(SDG805);
SDG805.InputBufferSize = 65535;
SDG805.OutputBufferSize = 65535;
if(~strcmp(SDG805.Status,'open'))
    warning('SDG805 port NOT open! Attempting....');
    fopen(SDG805);
end

switch(lower(mode))
    case 'on'
        modestr = 'ON';
    case 'off'
        modestr = 'OFF';
    otherwise
        error('Unknown mode string. Either ON or OFF!!')
end

scpi_sequence(SDG805, ...
    ['C1:OUTP ' modestr]);

end

