function SDG805_Pulsed(SDG805,PRF,PW,FRQ)
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

PRI = 1/PRF
per_frq = 1/FRQ;

ncyc = floor(PW/per_frq);



if(PRF > 0)
    %SDG805_OnOff(SDG805,'off');
    scpi_sequence(SDG805, ...
        'C1:BTWV STATE,ON', ...
        ['C1:BTWV PRD,' num2str(PRI)], ...
        ['C1:BTWV GATE_NCYC,NCYC'], ...
        ['C1:BTWV TIME,' num2str(ncyc)],'C1:BTWV STATE,ON');
    
    %SDG805_OnOff(SDG805,'on');
else
    
    scpi_sequence(SDG805, ...
        'C1:BTWV STATE,OFF');
end

end