function y = convert_16_to_2x8(u,endian)
% Function designed to convert a vector of 16bit values to a vector of 8bit
% values (by byte)

Nin = length(u);
Nout = 2*Nin;

y = zeros(1,Nout,'uint8');

% do MSB first
xMSB = bitshift(u,-8);
xLSB = bitand(u,255);

switch lower(endian)
    case 'big'
        y(1:2:end) = xMSB;
        y(2:2:end) = xLSB;
    case 'little'
        y(1:2:end) = xLSB;
        y(2:2:end) = xMSB;
    otherwise 
        error(['Invalid endian passed: ' endian '. Valid options are big or little']);
end

end