function [poles,zers] = polezero(num,den,k)
% Plot (or output) poles and zeros for a continuous transfer function
% Y(s)   K * (num)
% ---- = ------------
% U(s)     den
%
% Syntax: [poles,zers] = polezero(num,den,k)
%          (optional)
%         num  = numerator coefficients (s-domain)
%         den  = denominator coefficients (s-domain)
%         k    = gain (use 1 if embedded in num)
% Note: function will plot results if no output is requested
%
% K. Sawmiller, 2011

if(nargin == 2)
    k = 1.0;
end

poles = roots(den);
zers  = roots(k*num);

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
    %grid on;
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    title('Pole Zero Map');
    % axis limits
    xlim(1.1*[-max(abs(xlim)) max(abs(xlim))]);
    ylim(1.1*[-max(abs(ylim)) max(abs(ylim))]);
    xl = xlim(gca);
    yl = ylim(gca);
    hline = line([xl(1) xl(2)],[0 0]);
    hline = [hline line([0 0],[yl(1) yl(2)])];
    set(hline,'Color','k','LineStyle',':');
end

end