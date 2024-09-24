close all; clear all; clc

Fs = 10e6;
N = 60;
BW = 0.4e6;

for Fc = linspace(1.2e6,4.5e6,30)

h_BPF = fir_bpf_rect(Fc,BW,Fs,N);
h_plot = h_BPF;

if(~exist('hfig'))
    hfig = figure;
    subplot(1,2,1);
    th = linspace(0,2*pi,1000);
    plot(cos(th),sin(th),'k--');
    hold on; axis equal
    hp = plot(real(roots(h_plot)),imag(roots(h_plot)),'o');
    
    subplot(1,2,2)
    NFFT = 1000;
    df = Fs/NFFT;
    f = [0:df:(Fs-df)] - Fs/2;
    hp2 = plot(f, 10*log10(abs(fftshift(fft(h_plot,NFFT)))));
    title(['F_c ',num2str(Fc),'  ',num2str(N)])
else
    set(hp,'xdata',real(roots(h_plot)),'ydata',imag(roots(h_plot)))
    set(hp2,'ydata',10*log10(abs(fftshift(fft(h_plot,NFFT)))))
end
pause(0.2)
end