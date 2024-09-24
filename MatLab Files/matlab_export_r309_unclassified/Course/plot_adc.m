t = 0:0.01:10;
f = .2;
y = cos(2*pi*f*t);

fy = fft(y)./length(t);


delay = 1;

addpath('C:\Users\<user>\Documents\MATLAB\MATLAB\DSP');
addpath('C:\Users\<user>\Documents\MATLAB\MATLAB\Printing');
pp=PrepForPrint;
nBits = 1:14;
filename = 'adc_anim.gif';
filename2 = 'adc_anim2.gif';
for k = 1 : length(nBits)
    ys = ADC3(y,nBits(k));
    fys = fft(ys)./length(y);
    
    str = ['# of Bits: ' num2str(nBits(k))];
    if(k == 1)
        hf1 = figure;
        hp=plot(t,[y.' ys.']);
        ah = gca;
        grid off;
        ylim([-1.1 1.1]);
        ht = title(str);
        PrepForPrint(hf1,pp);
        set(hf1,'Color','w');
        drawnow;
        frame = getframe(hf1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        
        hf2 = figure;
        hp2=plot(t,10*log10(abs([fy.' fys.'])));
        ah2 = gca;
        grid off;
        ht2 = title(str);
        PrepForPrint(hf2,pp);
        set(hf2,'Color','w');
        xlim([0 0.5]);
        ylabel('Logarithmic Scale');
        xlabel('Frequency');
        drawnow;
        frame = getframe(hf2);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
        
 
        
    else
        set(hp(2),'YData',ys);
        set(ht,'String',str);
        drawnow;
        frame = getframe(hf1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay);
        
        set(hp2(2),'YData',10*log10(abs(fys)));
        set(ht2,'String',str);
        drawnow;
        frame = getframe(hf2);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename2,'gif','WriteMode','append','DelayTime',delay);
  
        
    end
    
    
    
    
    
end



