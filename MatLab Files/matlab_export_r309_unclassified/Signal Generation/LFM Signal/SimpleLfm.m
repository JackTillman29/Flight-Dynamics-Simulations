function [y,t,freq] = SimpleLfm(fStart,fStop,per,dt)
% help: [y,t,freq] = SimpleLfm(fStart,fStop,duration,dt)
    t = 0:dt:per-dt;
    wStart = 2*pi*fStart;
    wStop  = 2*pi*fStop;
    wRate  = (wStop - wStart) / per;
    freq = wStart + wRate.*t;
    y = exp(1j*(wStart.*t + 0.5*wRate.*t.^2));
end
