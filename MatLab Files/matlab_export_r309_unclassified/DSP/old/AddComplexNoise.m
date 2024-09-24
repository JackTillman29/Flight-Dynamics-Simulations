function y = AddComplexNoise(u,Vp)
    Vrms = single(sqrt(2)./2 * Vp);
    phase = single(2*pi*rand(1,length(u)));
    y = Vrms * exp(1i*phase) + u;
end