function p = PSKPulse( f, sbw, nBits, ovsmp, code)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% 13bit barker = [1 1 1 1 1 0 0 1 1 0 1 0 1]

if ( nargin == 3 )
    ovsmp = 10;
end

if ( nargin == 4 )
    code = round(rand(1,nBits));
end

subPulse = SquarePulse(f,sbw,ovsmp);

shifter = exp(1i*pi);

% build the command
s = 'out = [ ';
for k = 1:nBits
    if ( code(k) == 0 )
        s = [s 'subPulse '];
    elseif ( code(k) == 1 )
        s = [s 'subPulse * shifter '];
    end
end
s = [s '];'];
eval(s);
p = out;
end

