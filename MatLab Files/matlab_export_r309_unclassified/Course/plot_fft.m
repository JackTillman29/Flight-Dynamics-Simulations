Fs = 1000;
t = 0:(1/Fs):1;
df = Fs / length(t);
%y = cos(2*pi*5*t) + 0.1*cos(2*pi*25*t);
y = 0*t;
y(100:120) = 1;
y(200:240) = 0.8;
y(300:400) = 0.6;
y(500:700) = 0.5;
fy = fft(y)./length(t);
frqy = linspace(-Fs/2,Fs/2,length(t));
frqy = (0:length(t)-1) * Fs/length(t);


filename = 'fft.gif';

close all;
vidobj = VideoWriter('fft.avi','MPEG-4');
open(vidobj);

term = 0*t;
for k = 1 : ((length(frqy)-1)/2+1)
    if(k == 1)
        term = 0*t+fy(1);
    else
        k1 = k;
        k2 = length(fy)-k+2;
        term = term + ...
            exp(1i*2*pi*frqy(k1)*t) .* fy(k1) + ...
            exp(1i*2*pi*frqy(k2)*t) .* fy(k2);
    end
    if(k == 1)
        figure;
        hp=plot(t,real(term),'Color','r','LineWidth',2);
        hold on;
        plot(t,real(y),'Color',0.8*[1 1 1],'Linewidth',1);
        ylim([-0.2 1.2]);
        ht=title('1 FFT Term');
        set(gcf,'Color','w');
        set(gca,'XTickLabel','','YTickLabel','');
    else
        set(hp,'YData',real(term));
    end
    set(ht,'String',[num2str(k) ' FFT Terms']);
    drawnow;
    if(k == 1)
    pause;
    end
    currFrame = getframe(gcf);
    writeVideo(vidobj,currFrame);
        
    
end
close(vidobj);