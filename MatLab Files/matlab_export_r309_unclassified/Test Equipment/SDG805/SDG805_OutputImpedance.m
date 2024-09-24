function SDG805_OutputImpedance(obj,str)
switch(lower(str))
    case 'high-z'
        scpi_sequence(obj,'C1:OUTP LOAD,HZ');
    case '50'
        scpi_sequence(obj,'C1:OUTP LOAD,50');
    otherwise
        error(['Unknown output impedance!! ' str]);
end
end