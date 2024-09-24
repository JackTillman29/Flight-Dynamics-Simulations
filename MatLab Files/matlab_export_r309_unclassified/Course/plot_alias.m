Fs = 1000;
t = 0:(1/Fs):1;
decfac = 100;
Fs2 = Fs / decfac;
t2 = t(1:decfac:end);
delay = 0;
fArray = linspace(0,20,400);

filename = 'aliasing.mp4';

close all;

for k = 1 : length(fArray)

f = fArray(k);
f2 = mod(f,Fs/decfac);
y = cos(2*pi*f*t);
y2 = y(1:decfac:end);
y2_smooth = cos(2*pi*f2*t);

if(k == 1)
    hf=figure;
    subplot(1,2,1);
    hL1=plot(t,y);hold on;
    set(hL1,'Color',0.8*[1 1 1],'LineWidth',3);
    hLS = plot(t,y2_smooth);
    hL2=plot(t2,y2);hold off;
    set(hL2,'LineWidth',2);
    set(hL2,'Marker','o','MarkerEdgeColor','k', ...
        'MarkerSize',8,'MarkerFaceColor','r','LineStyle','none');
    ylim([-1.1 1.1]);
    set(hLS,'LineWidth',2);
    title('Discrete Time Domain');
    subplot(1,2,2);
    hLF1=plot(-5:5,fftshift(abs(fft(y2))/length(y2)));
    ylim([-1.1 1.1]);
    
    set(hLF1,'Marker','o','MarkerEdgeColor','k', ...
        'MarkerSize',8,'MarkerFaceColor','r','LineStyle','-');
    set(hLF1,'LineWidth',2);
    set(hLF1,'Color',[0.8500    0.3250    0.0980]);
    set(gca,'XTick',[-5 0 5],'XTickLabel',{'-F_s/2','0','F_s/2'});
    set(gcf,'Color','w');
    title('Discrete Frequency Domain');
    vidobj = VideoWriter(filename,'MPEG-4');
    vidobj.Quality = 15;
    vidobj.FrameRate = 10;
    open(vidobj);
     pause;
     
    
else
    set(hL1,'YData',y);
    set(hL2,'YData',y2);
    set(hLS,'YData',y2_smooth);
    
    set(hLF1,'YData',fftshift(abs(fft(y2))./length(y2)));
end
drawnow;
vidobj.writeVideo(getframe(gcf));

end
close(vidobj);