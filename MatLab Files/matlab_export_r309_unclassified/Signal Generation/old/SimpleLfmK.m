function [y,k] = SimpleLfmK(fracStart,fracStop,N)
% help: function [y,t,freq] = SimpleLfmK(fracStart,fracStop,N)
    k = 0:(N-1);
    wStart = 2*pi*fracStart;
    wStop  = 2*pi*fracStop;
    wRate  = (wStop - wStart)./N;
    %freq = wStart + wRate.*k;
    y = exp(1j*(wStart.*k + 0.5*wRate.*k.^2));
    
    
end