function [pw,prf,prf2,pwr] = GetWaveformProps(pulserise,hardzero,stream,dt)

    

    stream(stream < hardzero) = 0;
    stream(stream <= pulserise) = 0;
    
    stream_ref = stream;
    
    stream(stream ~= 0) = 1;
    
    rise_idx = find(diff(stream)==1);
    %fall_idx = find(diff(stream)==-1);

    pri = mean(diff(rise_idx)) * dt;
    prf = 1.0 / pri;
    
    idx.firstzero = find(stream == 0,1,'first');
    idx.firstnonzero = find(stream(idx.firstzero:end) > 0,1,'first') + idx.firstzero - 1;
    idx.secondzero = find(stream(idx.firstnonzero:end) == 0,1,'first') + idx.firstnonzero - 1;
    idx.secondnonzero = find(stream(idx.secondzero:end) > 0,1,'first') + idx.secondzero - 1;
    idx.thirdzero = find(stream(idx.secondnonzero:end) == 0,1,'first') + idx.secondnonzero - 1;
    idx.thirdnonzero = find(stream(idx.thirdzero:end) > 0,1,'first') + idx.thirdzero - 1;
    
    span = ((idx.secondzero-1) - idx.firstnonzero);
    pw = dt * span;
    pwr = sqrt(mean(stream_ref((idx.firstnonzero):(idx.secondzero))));
    
    if(isempty(prf))
        prf = 0;
    end
    fstream = abs(fft(single(stream)-mean(stream)));
    maxi = find(fstream > 1e11,1,'first');
    if(isempty(maxi))
        prf2 = 0;
    else
        prf2 = 1/dt * maxi / length(stream);
    end
    
end