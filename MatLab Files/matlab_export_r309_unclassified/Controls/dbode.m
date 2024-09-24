function [output,Fd] = dbode(numd,dend,kd,Ts,df,s2zmethod)
% Plot (or output) bode plot data for a discrete transfer function
% Yk   Kd * (numd)
% -- = ------------
% Uk     dend
%
% Syntax: [output,Fd] = dbode(numd,dend,kd,Ts,df)
%          (optional)
%         numd = numerator coefficients (z-domain)
%         dend = denominator coefficients (z-domain)
%         kd   = gain (use 1 if embedded in numd)
%         Ts   = Discrete sample time
%         df   = Frequency points to evaluate between 0:Nyquest
% Note: function will plot results if no output is requested
%
% K. Sawmiller, 2011

if(~exist('s2zmethod','var'))
    s2zmethod = 'tustin';
end

nyquistFreq  = 0.5 * 1.0 / Ts;  % Maximum for discrete bode
omega_vector = 1i * 2.0 * pi * (0:df:nyquistFreq);
omega_vector_z = s2z(omega_vector,Ts);
Fd = imag(omega_vector)./(2.0 * pi);
output = tfeval(numd,dend,kd,omega_vector_z);

if(nargout == 0)
    figure;
    subplot(2,1,1);
    hline = semilogx(Fd,20*log10(abs(output)));
    set(hline,'LineWidth',2.0);
    grid on;
    xlabel('Input Frequency (Hz)');
    ylabel('Magnitude (dB, 20 Log10)');
    title('Discrete Bode Plot, Magnitude');
    %axis tight;
    
    subplot(2,1,2);
    hline = semilogx(Fd,180/pi * unwrap(angle(output)));
    set(hline,'LineWidth',2.0);
    grid on;
    xlabel('Input Frequency (Hz)');
    ylabel('Output/Input Phase (deg)');
    title('Discrete Bode Plot, Phase');
    %axis tight;
end

end