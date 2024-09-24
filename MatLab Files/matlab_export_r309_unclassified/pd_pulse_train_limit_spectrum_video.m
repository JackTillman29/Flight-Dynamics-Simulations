close all; clear all; clc

addpath('Z:\sawmillkd\MATLAB\Printing\')
addpath('Z:\sawmillkd\MATLAB\PulsesAndWaveforms\')

pp = PrepForPrint;
pp.lineWidth = 1;

Fs = 30e6;
dt = 1/Fs;



IF  = 3e6;
PRI = 25e-6;
PW  = PRI/10;
N   = 100;


t = [0:dt:N*PRI];


L = length(t);


% carrier frequency
xc = exp(1i*2*pi*IF.*t);
% baseband pulse train
xb = PulseTrain(1/PRI,PW,L,dt);


% IF pulse train
x = xb.*xc;




figure('Position',[1217         529         560         420])
subplot(2,1,1)
plot(t,real(x))
xlim([0 3*PRI])


fx = fft(x);
mfx = abs(fx);
subplot(2,1,2)
plot(mfx)
% hold on
% plot(idx_max,mfx(idx_max),'r.')


idx_max = find(mfx == max(mfx),1);






%%
start_spread = 2000;
final_spread = 10;
spread_step = 1;
spread = [start_spread:(-1*spread_step):final_spread];
idx = 1:length(fx);

fund_spr = 9;

% spread = [2000 1000 500 100 10]



numPulses = 2;
% figure('Position',[1          62        1920         946])
% figure('Position',[1          62        600         500])
figure('Position',[10    62   981   688])

add_print_callbacks
% subplot(2,1,1)
% plot(t,real(x))
% xlim([0 numPulses*PRI])
% 
% fx = fft(x);
% mfx = abs(fx);
% subplot(2,1,2)
% plot(mfx)
% hold on
% plot(idx_max,mfx(idx_max),'r.')


make_movie = 1;
if(make_movie)
% aviobj = avifile('pd_windowed_spectrum_redux.avi','fps',15);
filename = 'pd_windowed_spectrum_redux_wo_fundamental.avi';
vidobj = VideoWriter(filename);
vidobj.FrameRate = 60;
open(vidobj)
end


iframe = 0;
for i = spread


    clc
    iframe = iframe + 1;
    disp(['frame: ',num2str(iframe),' of ',num2str(length(spread))])
    % compute FFT
    fx = fft(x);
    mfx = abs(fx);
    
    idx_max = find(mfx == max(mfx),1);
    sel_idx   = [ (idx_max-i) : (idx_max+i) ];
    zero_idx  = ~ismember(idx,sel_idx);
    zero_idx  = idx.*zero_idx;
    zero_idx(find(zero_idx == 0)) = [];
    
    % zero frequency components outside the frequency window
    fx(zero_idx) = 0;
    fx((idx_max-fund_spr):(idx_max+fund_spr)) = 0;
    mfx = abs(fx);
    
    xx = ifft(fx);
    
    subplot(2,1,1)
    plot(t,real(xx))
%     plot(t,10*log10(xx.^2))
    xlim([0 numPulses*PRI])
%     ylim([-60 10])
    ylim([-2 2])
%     title(['window length:',num2str(i*2)])
    title({'PD PRF Lines Filtered';...
        ['window length (centered around fundamental):',num2str(i*2)]})
    xlabel('Time [sec]')
    ylabel('Amplitude [V]')
    
    subplot(2,1,2)
%     plot(mfx)
    plot(10*log10(mfx))
    xlim([idx_max-start_spread idx_max+start_spread])
    ylim([-10 40])
    xlabel('Freq Bin Number')
    ylabel('Counts')
    
    if(make_movie)
        writeVideo(vidobj,getframe(gcf))
    else
        pause(0.01)
        PrepForPrint(gcf,pp)
%         filename = ['pd_filter_prf_lines_',num2str(i)];
%         saveas(gcf,[filename '.fig'])
%         kprint(get(gcf,'Number'),[filename '.png'],180)
    end
    
    
    
    
    
    






end


if(make_movie)
close(vidobj)
end