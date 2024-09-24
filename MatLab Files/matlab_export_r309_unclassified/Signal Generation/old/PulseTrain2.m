function [outdat] = PulseTrain2(varargin)
% Help: 
%         'PRF'             Hz
%         'PulseDuration'   sec
%         'TrainDuration'   sec (use this or NumberOfPulses or MaxSamples)
%         'NumberOfPulses'  #
%         'SampleTime'      sec
%         'SampleRate'      Hz (use this or SampleTime)
%         'MaxSamples'      # length of output (will truncate)
%         'Carrier'         Hz carrier applied (if nonzero)
%         'Modulation'      Hz modulation applied (if nonzero)
%         'Delay'           sec delay from zero seconds
%         'Type'            'Single' if single precision is desired


input.PRF = 1;
input.PW  = 0.1;
input.Np  = 10;
input.dt  = 0.05;
input.LenLim = 0;
input.Carrier = 0;
input.Delay = 0;
input.TrainDuration = 0;
input.Type = 'double';

for k = 1:2:length(varargin)
    thisString = varargin{k};
    thisValue  = varargin{k+1};
    switch thisString
        case 'PRF'
            input.PRF = thisValue;
        case 'PulseDuration'
            input.PW = thisValue;
        case 'TrainDuration'
            input.TrainDuration = thisValue;
        case 'NumberOfPulses'
            input.Np = thisValue;
        case 'SampleTime'
            input.dt = thisValue;
        case 'SampleRate'
            input.dt = 1.0 / thisValue;
        case 'MaxSamples'
            input.LenLim = thisValue;
        case 'Carrier'
            input.Carrier = thisValue;
        case 'Modulation'
            input.Modulation = thisValue;
        case 'Delay'
            input.Delay = thisValue;
        case 'Type'
            input.Type = thisValue;
    end
end

% compute reference data
pri = 1.0 / input.PRF;

%error check
if(pri < input.dt)
    error('PRI of waveform is LESS than DT!');
end

if(input.TrainDuration ~= 0)
    lenSamples = round(input.TrainDuration / input.dt);
    input.Np = round(input.TrainDuration / pri);
else
    lenSamples = round(input.Np * pri / input.dt);
end
lenPri     = round(pri / input.dt);
lenPulse   = round(input.PW / input.dt);

tSamples = (1:lenSamples) - 1;

outdat.waveform = (mod(tSamples,lenPri) <= lenPulse);
outdat.time = input.dt * tSamples;
outdat.carrier = exp(1i*2*pi*input.Carrier*outdat.time);
if(input.LenLim)
    outdat.waveform = outdat.waveform(1:input.LenLim);
    outdat.time = outdat.time(1:input.LenLim);
    outdat.carrier = outdat.carrier(1:input.LenLim);
end

if(input.Carrier ~= 0)
    if((1.0 / input.Carrier) < input.dt)
        error(['Carrier too high! Max supported at this Fs is ' num2str(1.0/input.dt)]);
    end
    outdat.waveform = outdat.waveform .* exp(1i*2*pi*input.Carrier*outdat.time);
end

if(input.Modulation ~= 0)
    outdat.waveform = outdat.waveform .* exp(1i*2*pi*input.Modulation*outdat.time);
end

if(input.Delay ~= 0)
    sample_delays = round(input.Delay / input.dt);
    outdat.waveform((sample_delays+1):end) = outdat.waveform(1:(end-sample_delays));
    outdat.waveform(1:sample_delays) = 0;
end

if(nargout == 0)
    figure;
    plot(outdat.time,([real(outdat.waveform)' abs(outdat.waveform)' real(outdat.carrier)']));
    ylim([-1 2]);
    title('Pulsetrain');
    legend('R(sig)','Env(sig)','R(carrier)');
end

if(strcmpi(input.Type,'Single'))
    outdat.time = single(outdat.time);
    outdat.waveform = single(outdat.waveform);
    outdat.carrier = single(outdat.carrier);
end

end