function y = ApplyPhaseModulation(u,Ts,code_in,chipLen_s)
% us is input signal (typically unmodulated)
% code is phase code in radians, or a string of 1's and 0's
% chipLen is length of chip in seconds

% first, determine the number of discrete samples for the phase shifting
chipLen_N = round(chipLen_s / Ts);

% if a string was passed, convert it to radians +/- pi/2
if(ischar(code_in))
    code = zeros(1,length(code_in));
    for k = 1 : length(code_in)
        code(k) = (2*str2double(code_in(k))-1) .* (pi/2);
    end
else
    code = code_in;
end

y = u;
idx1 = 1;

% for all chips
for k = 1 : length(code)
    idx2 = idx1 - 1 + chipLen_N;
    if(idx2 > length(u))
        idx2 = length(u);
        warning(['Truncating phase code on chip # ' num2str(k)]);
    end
    y(idx1:idx2) = y(idx1:idx2) .* exp(1i*code(k));
    idx1 = idx2 + 1;
end
end