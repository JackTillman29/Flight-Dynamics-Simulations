for k = 1 : length(tgtRanges)
    tgtRange = tgtRanges(k);
    tgtAmplitude = tgtAmplitudes(k);
    tgtDoppler = tgtDopplers(k);

    % generate impulse return for this target
    rxImpulse = gen_pd_impulse(...
        rxImpulse,tgtRange,PRF,...
        Ts,tgtAmplitude,rxLO,...
        rxRF,tgtDoppler);
end

rxSamples = conv(rxImpulse,txWaveform,'same');
