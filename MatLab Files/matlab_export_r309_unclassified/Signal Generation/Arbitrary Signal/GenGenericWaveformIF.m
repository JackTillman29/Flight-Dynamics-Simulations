function wav = GenGenericWaveform(...
    timeTableF, freqTable, timeTableA, ampTable, Ts, varargin)
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
% VARARGIN:
%   the first arg can be a number, N
%      N is an optional argument that will force the output to be length N
%   otherwise, the rest of the inputs follow the String-Value pair
%   convention. Example:
%        <'Color',[0 0 1]>
%    e.g. plot(x,y,'Color',[0 0 1])
%
% String-Value pairs:
%   'rf' and 'bw' inputs will automatically determine which parts of a
%   wideband signal are observable (assumes perfect rectangular IF filter*)
%       'rf' = defines the center RF frequency of the band-of-interest
%       'bw' = defines the bandwidth of the band-of-interest
%
%   'if_filter_freq' and 'if_filter_gain' defines the frequency / gain
%   character of the IF filter.
%       'if_filter_freq' = array of baseband freq points (centered over
%           band-of-interest in the code)
%           ****in Hz
%       'if_filter_gain' = array of gains corresponding to each freq point
%           in 'if_filter_freq'
%           ****in MAGNITUDE. [not in dB]
%
%  TODO: add 'if_filter_gain' and 'if_filter_phase' parameters so that a
%       lookup will occur and modify the original interpolated ampTable
%       values. (DONE, now need to add phs lookup)

% parse varargin
if(length(varargin)>0)
    parsingArgs = 1;
    k = 1;
    while(parsingArgs)
        if(isnumeric(varargin{k}))
            N = varargin{k};
        elseif(ischar(varargin{k}))
            % string-value pair
            prop = varargin{k};
            k = k + 1;
            val  = varargin{k};
            switch lower(prop)
                case 'rf'
                    rf = val;
                case 'bw'
                    bw = val;
                case 'if_filter_freq'
                    if_filter_freq = val;
                case 'if_filter_gain'
                    if_filter_gain = val;
                case 'if_filter_phs'
                    if_filter_phs = val;
                otherwise
                    error('unknown property')
            end
        end
        k = k + 1;
        % done?
        if(k > length(varargin))
            parsingArgs = 0;
        end
    end
end

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

% Extract band-of-interest and modify interpolated freq / amp arrays
if(exist('rf') & exist('bw'))
    fLower = rf - bw/2; fUpper = rf + bw/2;
    idx = find( (fInterp >= fLower) & (fInterp <= fUpper) );
    aInterpNew = zeros(size(aInterp));
    if(~exist('if_filter_gain'))
        % perfect rectangular IF filter
        aInterpNew(idx) = aInterp(idx);
    else
        % TODO: include 'if_filter_phs' term (constants added after phase
        %                                   integration)
        % set to zero so interp1 extrapolates to zero outside bounds of
        % input freq array
        stopbandGain    = 0;
        aInterpNew(idx) = aInterp(idx) .* interp1(...
            rf+if_filter_freq, if_filter_gain,fInterp(idx),...
            'linear',stopbandGain);
    end
    aInterp = aInterpNew;
end

% Integrate over frequency function to generate waveform
wav = aInterp .* exp(1i*cumsum(2*pi*fInterp).*Ts);

end