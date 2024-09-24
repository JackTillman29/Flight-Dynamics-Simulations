function [yout] = step(num,den,t)
error('not quite working yet. Use simulink or fix this and send me the updates!');
% Plot (or output) step response for a continuous transfer function
% Y(s)   K * (num)
% ---- = ------------
% U(s)     den
%
% Syntax: [uout,yout] = step(num,den,t)
%          (optional)
%         num  = numerator coefficients (s-domain)
%         den  = denominator coefficients (s-domain)
%         t    = time vector 0:0.1:10
%
%
% K. Sawmiller, 2011
    n = length(t);
    yout = zeros(1,n);
    
    % Convert transfer function to step response (*1/s)
    tf(1).num = num;
    tf(1).den = den;
    tf(2).num =  1;
    tf(2).den = [1 0];
    
    % Linearly combine the transfer functions into a single composite tf
    combtf = tf_combine(tf);
    
    % Decompose the transfer function into a summation of partial fractions
    [partNums,partDens,partConsts] = residue(combtf.num,combtf.den);
    disp(partNums);
    disp(partDens);
    disp(partConsts);
    
    % Now we have to loop through the partial fractions to expand into a
    % more preferrable form. MATLAB handles repeated roots funny.
    p = 1;
    while(p <= length(partNums))
        m = find(partNums(p:end)~=0,1,'first') - 1;
        p = p + m;
        disp(['Multiplicity: ' num2str(m)]);
        disp(['This P = ' num2str(p)'])
    end
    
    
    if(nargout == 0)
        figure;
        plot(t,yout);
        xlabel('Time (sec)');
        ylabel('Response');
        title('Continuous Transfer Function Step Response');
        legend('Input','Output');
        grid on;
    end
end