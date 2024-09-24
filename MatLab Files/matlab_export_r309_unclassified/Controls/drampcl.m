function [uout,yout,rout] = drampcl(num,den,n,start,slope,respmafn)
% Plot (or output) ramp response for a discrete transfer function
% Y(z)   Kd * (numd)
% ---- = ------------
% U(z)     dend
%
% Syntax: [uout,yout] = stepd(num,den,n,start)
%          (optional)
%         numd  = numerator coefficients (z-domain)
%         dend  = denominator coefficients (z-domain)
%         n     = number of frames to simulate
%         start = frame to start ramp input
%         slope = ramp slope
%
%   NOTE: Requires denominator leading term be 1
%
% K. Sawmiller, 2011

maf = zeros(1,respmafn);
imaf = 1;

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
    rout = zeros(1,n);
    qout = uout;
    
    for k = 1:n
        rout(k) = mean(maf);
        if(k > start)
            qout(k) = slope*(k+1-start);
            u(1) = qout(k) - rout(k);
            %u(1) = 1.0;
        else
            u(1) = 0.0;
        end
        
        y(1) = sum(num.*u) - sum(den(2:end).*y(2:end));
        maf(imaf) = y(1);
        disp(imaf)
        imaf = imaf+1;
        if(imaf > respmafn)
            imaf = 1;
        end
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
        plot([rout' qout'],'.-');
        xlabel('Sample #');
        ylabel('Magnitude (Linear)');
        title('Discrete Transfer Function Step Response');
        legend('Output');
        grid on;
    
end