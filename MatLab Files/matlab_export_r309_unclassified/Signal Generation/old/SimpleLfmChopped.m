function [y,t,freq] = SimpleLfmChopped(fStart,fStop,per,dt,fc,narrowBW)
% help: [y,t,freq] = SimpleLfm(fStart,fStop,duration,dt)
    t = 0:dt:per;
    wStart = 2*pi*fStart;
    wStop  = 2*pi*fStop;
    wRate  = (wStop - wStart) / per;
    freq = wStart + wRate.*t;
    freqHz = freq/(2*pi); 
    centerOfLFM = intersect(find(freqHz>=fc-narrowBW*.5),...
        find(freqHz<=fc+narrowBW*.5) );

    y = exp(1j*(wStart.*t(centerOfLFM) ...
        + 0.5*wRate.*t(centerOfLFM).^2));
    
    t = 0:dt:length(centerOfLFM)*dt; 
end
