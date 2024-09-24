close all;
addpath('I:\MATLAB\Printing');
addpath('I:\MATLAB\Windows');

N = 1000;
Nk = 10;
u = randn(1,N)+randn(1,N)*1i;
kernel = zeros(1,N);
kernel(1:10) = 1;

finput = fft(u)./N;
fkernel = fft(kernel)./Nk;
foutput = finput .* fkernel;
y = ifft(foutput);

% Create more agressive FIR filter
fvec = linspace(-.5,0.5,N);
firWindow = double(abs(fvec) <= 0.1);
firWindow(firWindow == 0) = 1e-6;
firWindow = fftshift(firWindow);

foutput_fir = finput .* firWindow;


figure;
plot(fvec,20*log10(abs(fftshift([fkernel.' finput.' firWindow']))));
ylim([-60 10]);
xlabel('x Fs');
ylabel('20log_{10}');
title('Frequency Domain - PreConvolution');
grid on;legend('FIR','Pulse','Input');
add_print_callbacks;


figure;
plot(fvec,20*log10(abs(fftshift([foutput_fir.' foutput.' finput.' ]))));
ylim([-60 10]);
xlabel('x Fs');
ylabel('20log_{10}');
title('Frequency Domain - PostConvolution');
grid on;legend('Input','FIR','Pulse');
add_print_callbacks;

cpsd_input = cumsum(fftshift(finput .* conj(finput)));
cpsd_output = cumsum(fftshift(foutput .* conj(foutput)));
cpsd_output_fir = cumsum(fftshift(foutput_fir .* conj(foutput_fir)));

figure;
plot(fvec,[cpsd_output_fir.' cpsd_output.' cpsd_input.']./cpsd_input(end));
xlabel('x Fs');
ylabel('CPSD Normalized (|x|^2)')
grid on;title('Cumulative Power Spectra');
legend('FIR','Pulse','Input');
ylim([-.1 1.1]);
add_print_callbacks;

%% what happens with a pulse in the return???
signal_return = u+10*fftshift(kernel); % simply drop a signal in the return PRI
fsignal = fft(signal_return)./N;

fsignal_output = fsignal .* fkernel;
fsignal_output_fir = fsignal .* firWindow;

figure;
hold on;
plot(fvec,20*log10(abs(fftshift(fsignal_output))));
plot(fvec,20*log10(abs(fftshift(fsignal_output_fir))));
xlabel('x Fs');
ylabel('20log_{10}')
grid on;title('Power Spectrum');
legend('Pulse','FIR');
add_print_callbacks;

figure;
hold on;
plot(20*log10(abs(N*ifft(fsignal_output))));
plot(20*log10(abs(N*ifft(fsignal_output_fir))));
plot(20*log10(abs(signal_return)));
xlabel('Sample');
ylabel('Amplitude')
grid on;title('Time Domain');
legend('Pulse','FIR','Input');
add_print_callbacks;


%% what about windowing the transmit kernel to reduce sidelobes??
kernel_wnd = zeros(1,N);
Nk = 17;
wnd_factor = 10.^(6/20);
kernel_wnd(1:Nk) = 1*wnd_factor;
kernel_wnd(1:Nk) = hamming(Nk).*kernel_wnd(1:Nk);
fkernel_wnd = fft(kernel_wnd)./Nk;

figure;
plot(fftshift([kernel.' kernel_wnd.']))

figure;
hold on;
plot(fvec,20*log10(abs(fftshift(fkernel))));
plot(fvec,20*log10(abs(fftshift(fkernel_wnd))));
ylim([-50 10]);
xlabel('x Fs');
ylabel('20log_{10}');
grid on;
title('Matched Filter');
legend('Transmit','Tuned Hamming');
add_print_callbacks;


%%
foutput_wnd = finput .* fkernel_wnd;

figure;
plot(fvec,20*log10(abs(fftshift([foutput_wnd.' foutput.' finput.' ]))));
ylim([-60 10]);
xlabel('x Fs');
ylabel('20log_{10}');
title('Frequency Domain - PostConvolution');
grid on;legend('Input','Hamming','Pulse');
add_print_callbacks;

cpsd_input = cumsum(fftshift(finput .* conj(finput)));
cpsd_output = cumsum(fftshift(foutput .* conj(foutput)));
cpsd_output_wnd = cumsum(fftshift(foutput_wnd .* conj(foutput_wnd)));

figure;
plot(fvec,[cpsd_output_wnd.' cpsd_output.' cpsd_input.']./cpsd_input(end));
xlabel('x Fs');
ylabel('CPSD Normalized (|x|^2)')
grid on;title('Cumulative Power Spectra');
legend('Hamming','Pulse','Input');
ylim([-.1 1.1]);
add_print_callbacks;


fsignal_output_wnd = fsignal .* fkernel_wnd;

figure;
hold on;
plot(fvec,20*log10(abs(fftshift(fsignal_output))));
plot(fvec,20*log10(abs(fftshift(fsignal_output_wnd))));
xlabel('x Fs');
ylabel('20log_{10}')
grid on;title('Power Spectrum');
legend('Pulse','Hamming');
add_print_callbacks;

figure;
hold on;
plot(20*log10(abs(N*ifft(fsignal_output))));
plot(20*log10(abs(N*ifft(fsignal_output_wnd))));
plot(20*log10(abs(signal_return)));
xlabel('Sample');
ylabel('Amplitude')
grid on;title('Time Domain');
legend('Pulse','Hamming','Input');
add_print_callbacks;