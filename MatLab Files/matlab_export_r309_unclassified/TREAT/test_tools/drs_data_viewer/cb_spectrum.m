function cb_spectrum( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

haxes = findobj('Type','Axes','Parent',gcf);

%for k = 1 : length(haxes)
    hax = haxes(1);
    xb=get(hax,'XLim');
    xvec = evalin('base','xvec;');
    dat1 = evalin('base','dat1;');
    dat2 = evalin('base','dat2;');
    iStart = find(xvec > xb(1),1,'first');
    iStop = find(xvec < xb(2),1,'last');
    
    ezja_header = evalin('base','ezja_header;');
    if(~ezja_header)
        i1_counts = dat1.pageMemory(iStart:iStop,1);
        q1_counts = dat1.pageMemory(iStart:iStop,2);
        i2_counts = dat2.pageMemory(iStart:iStop,1);
        q2_counts = dat2.pageMemory(iStart:iStop,2);

        v1 = single(i1_counts) + 1i*single(q1_counts);
        v2 = single(i2_counts) + 1i*single(q2_counts);
    else
         dat1.temp = single(dat1.pageMemory(iStart:iStop));
         dat2.temp = single(dat1.pageMemory(iStart:iStop));
         disp('Applying Hilbert Transform on Channel 1,2');
         dat1.tempf = fft(dat1.temp);
         dat2.tempf = fft(dat2.temp);
         dat1.tempf(floor(length(dat1.tempf)/2):end) = 0;
         dat2.tempf(floor(length(dat2.tempf)/2):end) = 0;
         v1 = ifft(dat1.tempf);
         v2 = ifft(dat2.tempf);
         disp('Done');
        
    end
    
    f1 = fftshift(fft(v1));
    f2 = fftshift(fft(v2));
    
    fvec = evalin('base','Fs;') * (0:length(f1)-1)./length(f1);
    fvec = fvec - fvec(end)/2;
    
    
    
    same = 0;
    if(same)
        figure;
        plot(1e-6*fvec,[abs(20*log10(abs(f1))) abs(20*log10(abs(f2)))]);
        grid on;
        xlabel('Frequency (MHz)');
        ylabel('Amplitude (20log_{10})');
        title('Spectrum');
        legend('CH1','CH2');
    else
        figure;
        subplot(2,1,1);
        plot(1e-6*fvec,[abs(20*log10(abs(f1)))]);
        grid on;
        xlabel('Frequency (MHz)');
        ylabel('Amplitude (20log_{10})');
        title('Spectrum [CH1]');ah = gca;
        
        subplot(2,1,2);
        plot(1e-6*fvec,[abs(20*log10(abs(f2)))]);
        grid on;
        xlabel('Frequency (MHz)');
        ylabel('Amplitude (20log_{10})');
        title('Spectrum [CH2]');ah = [ah gca];
        linkaxes(ah,'xy');
        
    end
    
%end


end

