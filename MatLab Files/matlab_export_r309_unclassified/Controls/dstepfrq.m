function [uout,yout] = dstepfrq(num,den,n,start,rad_per_sample) 
% Plot (or output) step response for a discrete transfer function
% Y(z)   Kd * (numd)
% ---- = ------------
% U(z)     dend
%
% Syntax: [uout,yout] = stepd(num,den,n,start)
%          (optional)
%         numd  = numerator coefficients (z-domain)
%         dend  = denominator coefficients (z-domain)
%         n     = number of frames to simulate
%         start = frame to start step input
%
%   NOTE: Requires denominator leading term be 1
%
% K. Sawmiller, 2011
lenDiscrep = length(den) - length(num);

if(lenDiscrep > 0) % denominator is higher order
    num = [zeros(1,lenDiscrep) num];
else
    den = [zeros(1,-lenDiscrep) den];
end
    
    u = zeros(1,length(num));
    y = zeros(1,length(den));
    
    % temporarily fix the dimensions of the num/den so they multiply ok
    
    
    uout = zeros(1,n);
    yout = zeros(1,n);
    
    for k = 1:n
        if(k > start)
            u(1) = 1.0 * exp(1i*k*rad_per_sample);
        else
            u(1) = 0.0;
        end
        
        y(1) = sum(num.*u) - sum(den(2:end).*y(2:end));
        yout(k) = y(1);
        uout(k) = u(1);
        
        % shift registers
        y(2:end) = y(1:(end-1));
        u(2:end) = u(1:(end-1));
        
    end
    disp(['Final u: ' num2str(u)]);
    disp(['Final y: ' num2str(y)]);
    
    if(nargout == 0)
        figure;
        plot([uout' yout']);
        xlabel('Sample #');
        ylabel('Magnitude (Linear)');
        title('Discrete Transfer Function Step Response');
        legend('Input','Output');
        grid on;
    
end