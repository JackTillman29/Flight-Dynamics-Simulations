function [Klead AmbigArray] = BPSKAmbiguity(Fc, Fs, PW, BitKey, Modulations)

dt = 1/Fs;
TimePulse = (0:dt:PW);
ThreatCarrier = SimpleLfm(Fc,Fc,PW,dt);
ThreatBPSKShifter = SimplePSK(BitKey,TimePulse);
ThreatPulse = ThreatCarrier.*exp(1i.*(ThreatBPSKShifter*pi));
RcvrMatch = ThreatPulse;

%Threat Waveform (Skin Return)
[match.ThreatPulse, klead.ThreatPulse] = PulseCompressionXCORRs(ThreatPulse, RcvrMatch);
AmbigArray = zeros(length(Modulations),length(klead.ThreatPulse));

for k = 1:length(Modulations)
    PulseMod = ThreatPulse.*exp(1i*2*pi*Modulations(k).*TimePulse);
    [match.PulseMod, klead.PulseMod] = PulseCompressionXCORRs(PulseMod, RcvrMatch);
    AmbigArray(k,:) = match.PulseMod;
end

Klead = klead.PulseMod;

end