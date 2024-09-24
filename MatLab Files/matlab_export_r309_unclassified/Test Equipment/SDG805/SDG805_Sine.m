function SDG805_Sine(SDG805,freq,magnitude,magtype)
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

if(~exist('magtype'))
    magtype = 'Vp';
    warning('Assuming input is in volts peak (amplitude)');
end

switch(lower(magtype))
    case 'vp'
        ampval = magnitude;
    case 'dbm'
        dbm = magnitude+6;
        PmW = 10^(dbm/10);
        vrms = sqrt(PmW * 50 / 1e3);
        ampval = sqrt(2) * vrms;
    otherwise
        error('Unknown magnitude type!!')
end


scpi_sequence(SDG805, ...
    'C1:BSWV WVTP,SINE', ...
    ['C1:BSWV FRQ,' num2str(freq)], ...
    ['C1:BSWV AMP,' num2str(ampval)] ...
    );

end