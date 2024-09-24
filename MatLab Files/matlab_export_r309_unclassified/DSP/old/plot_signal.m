function plot_signal(x,subtitle,hf,Fs,wnd,fScaleFact,fOffset,fLabel)

% calculate power loss from windowing and apply it back in
p_wnd = mean(wnd.^2); % this is the power factor associated with the window function

N = length(x);
df = Fs / N;

x_fft = fft(wnd.*x) ./ N;
x_fft_freq = Fs * ((0:(N-1)) ./ N);
x_psd = (abs(x_fft) .^ 2) / df / p_wnd;
x_cpsd = cumsum(x_psd) .* df;

figure(hf);
subplot(2,1,1);
plot( ...
    fScaleFact*(x_fft_freq-fOffset), 10*log10(x_psd) );
xlabel(fLabel);
ylabel('dBW/Hz');
title(['Power Spectral Density, ' subtitle]);
grid on;
ah = gca;


subplot(2,1,2);
plot( ...
    fScaleFact*(x_fft_freq-fOffset), x_cpsd ,'.-');
xlabel(fLabel);
ylabel('W');
title(['Cumulative Power Spectrum, ' subtitle]);
grid on;
ah = [ah gca];
linkaxes(ah,'x');

print_props = PrepForPrint();
print_props.xylabelFontSize = 12;
print_props.legendFontSize = 10;
PrepForPrint(hf,print_props);
end