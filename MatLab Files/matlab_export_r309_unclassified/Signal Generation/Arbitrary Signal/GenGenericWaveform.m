function wav = GenGenericWaveform(timeTableF,freqTable,timeTableA,ampTable,Ts,N)
% function wav = GenGenericWaveform(timeTable,freqTable,ampTable,Ts,<opt>N)
% K. Sawmiller 2018
% ex: Simple LFM waveform, 128us, 0 to +5MHz, amplitude 1, sampled at 1GHz
% GenGenericWaveform([0 128e-6], ... % Time Array (freq data)
%                    [0 5e6]   , ... % Frequency Array
%                    [0 128e-6], ... % Time Array (amp data)
%                    [1 1]     , ... % Amplitude Array
%                    1.0e-9    )    % Sample Time
%
% NOTE: The amplitude and frequency data are interpolated on their own
%       time arrays. This allows for complex amplitude & simple freq or
%       vice versa.
% 
% Same LFM waveform with sinusoidal AM @ 20kHz
% GenGenericWaveform([0 128e-6],
%       [0 5e6],
%       [0:.1:1]*128e-6,
%       sin(2e4*[0:.1:1]*128e-6),
%       Ts);
%
% N is an optional argument that will force the output to be length N

if( ...
        ~isa(timeTableF,'double') || ...
        ~isa(timeTableA,'double') || ...
        ~isa(Ts,'double') || ...
        ~isa(freqTable,'double') || ...
        ~isa(ampTable,'double') ...
        )
    error('routine does not work with single precision during generation. cast result as single if needed.');
end
% Create internal time reference vector based on table
maxTime = max(timeTableA(end),timeTableF(end));

if(~exist('N'))
    refTimeVec = (timeTableF(1):Ts:maxTime);
else
    refTimeVec = (0:(N-1)) * Ts;
end

% Compute interpolated instantaneous frequency
fInterp = interp1(timeTableF,freqTable,refTimeVec,'linear','extrap');

% Compute interpolated instantaneous amplitude
aInterp = interp1(timeTableA,ampTable,refTimeVec,'linear','extrap');

% Integrate over frequency function to generate waveform
wav = aInterp .* exp(1i*cumsum(2*pi*fInterp).*Ts);

end