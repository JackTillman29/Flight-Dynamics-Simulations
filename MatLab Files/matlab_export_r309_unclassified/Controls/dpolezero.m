function [poles,zers] = dpolezero(numd,dend,kd)
% Plot (or output) poles and zeros for a discrete transfer function
% Y(z)   Kd * (numd)
% ---- = ------------
% U(z)     dend
%
% Syntax: [poles,zers] = dpolezero(numd,dend,kd)
%          (optional)
%         numd  = numerator coefficients (z-domain)
%         dend  = denominator coefficients (z-domain)
%         kd    = gain (use 1 if embedded in num)% Note: function will plot results if no output is requested
%
% K. Sawmiller, 2011

if(nargin == 2)
    kd = 1.0;
end

poles = roots(dend);
zers  = roots(kd*numd);

if(nargout == 0)
    figure;
    h=plot(real(zers), imag(zers), 'go', ...
         real(poles),imag(poles),'bx');
    if(~isempty(zers) && ~isempty(poles))
        legend('Zeros','Poles');
    elseif(~isempty(poles))
        legend('Poles');
    else
        legend('Zeros');
    end
    set(h,'MarkerSize',12.0);
    % grid on;
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    title('Pole Zero Map');
    % axis limits
    maxval = max([abs(xlim) abs(ylim)]);
    xlim(1.1*[-max(abs(maxval)) max(abs(maxval))]);
    ylim(1.1*[-max(abs(maxval)) max(abs(maxval))]);
    %set(gca,'DataAspectRatio',[1 1 1]);
    theta = 0:0.001:(2*pi + 0.001);
    
    hline = line([cos(theta)],[sin(theta)]);
    set(hline,'Color','k','LineStyle',':');
    axis square;
end

end