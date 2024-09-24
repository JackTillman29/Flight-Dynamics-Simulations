close all; clear all; clc

Fs = 10e6;
N = 11;

% h_LPF = fir_lpf_rect(Fc,Fs,N,'plot');

for Fc = linspace(1e3,4e6,100)
    
h_LPF = fir_lpf_rect(Fc,Fs,N);
h_plot = h_LPF;

if(~exist('hfig'))
    hfig = figure;
    subplot(1,3,1);
    th = linspace(0,2*pi,1000);
    plot(cos(th),sin(th),'k--');
    hold on; axis equal
    hp = plot(real(roots(h_plot)),imag(roots(h_plot)),'o');
    xlim(5*[-1 1])
    ylim(5*[-1 1])
    
    subplot(1,3,2)
    NFFT = 1000;
    df = Fs/NFFT;
    f = [0:df:(Fs-df)] - Fs/2;
    hp2 = plot(f, 10*log10(abs(fftshift(fft(h_plot,NFFT)))));
    htitle = title(['F_c ',num2str(Fc),'  ',num2str(N)]);
    
    subplot(1,3,3)
    hp3 = stem(h_LPF);
else
    set(hp,'xdata',real(roots(h_plot)),'ydata',imag(roots(h_plot)))
    set(hp2,'ydata',10*log10(abs(fftshift(fft(h_plot,NFFT)))))
    set(hp3,'ydata',h_LPF)
    set(htitle,'String',['F_c ',num2str(Fc),'  ',num2str(N)])
end
pause(0.2)
end






