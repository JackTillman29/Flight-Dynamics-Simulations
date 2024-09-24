function [output,F] = bode(num,den,k,Fmax,df)
% Plot (or output) bode plot data for a discrete transfer function
% Y(s)   K * (num)
% ---- = ------------
% U(s)     den
%
% Syntax: [output,F] = bode(num,den,k,Fmax,df)
%          (optional)
%         num  = numerator coefficients (s-domain)
%         den  = denominator coefficients (s-domain)
%         k    = gain (use 1 if embedded in num)
%         Fmax = Maximum frequency to evaluate (Hz)
%         df   = Frequency points to evaluate between 0:Fmax
% Note: function will plot results if no output is requested
%
% K. Sawmiller, 2011

omega_vector = 1i * 2.0 * pi * (0:df:Fmax);
F = imag(omega_vector)./(2.0 * pi);
output = tfeval(num,den,k,omega_vector);

if(nargout == 0)
    figure;
    subplot(2,1,1);
    hline = semilogx(F,20*log10(abs(output)));
    set(hline,'LineWidth',2.0);
    grid on;
    xlabel('Input Frequency (Hz)');
    ylabel('Magnitude (dB, 20 Log10)');
    title('Bode Plot, Magnitude');
    %axis tight;
    
    subplot(2,1,2);
    hline = semilogx(F,180/pi * angle(output));
    set(hline,'LineWidth',2.0);
    grid on;
    xlabel('Input Frequency (Hz)');
    ylabel('Output/Input Phase (deg)');
    title('Bode Plot, Phase');
    %axis tight;
    
    figure;
    hline = plot(180/pi * unwrap(angle(output)), 20*log10(abs(output)));
    set(hline,'LineWidth',2.0);
    grid on;
    xlabel('Phase (deg)');
    ylabel('Magnitude (dB)');
    title('Nichols Plot');
    %axis tight;
end

end