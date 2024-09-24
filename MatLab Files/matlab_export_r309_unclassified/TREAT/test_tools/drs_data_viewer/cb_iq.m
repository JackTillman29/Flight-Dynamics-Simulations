function cb_iq( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

haxes = findobj('Type','Axes','Parent',gcf);
ezja_header = evalin('base','ezja_header;');

%for k = 1 : length(haxes)
    hax = haxes(1);
    xb=get(hax,'XLim');
    xvec = evalin('base','xvec;');
    dat1 = evalin('base','dat1;');
    dat2 = evalin('base','dat2;');
    iStart = find(xvec > xb(1),1,'first');
    iStop = find(xvec < xb(2),1,'last');
    
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
    
    
    
    ts = evalin('base','Ts;');
    

    figure;
    subplot(2,1,1);
    plot(ts*(0:length(v1)-1),[real(v1) imag(v1)],'.-');
    grid on;
    %xlabel('Frequency (MHz)');
    %ylabel('Amplitude (20log_{10})');
    title('I/Q [CH1]');ah = gca;
    xlabel('Time (sec)');

    subplot(2,1,2);
    plot(ts*(0:length(v1)-1),[real(v2) imag(v2)],'.-');
    grid on;
    %xlabel('Frequency (MHz)');
    %ylabel('Amplitude (20log_{10})');
    title('I/Q [CH2]');ah = [ah gca];

    xlabel('Time (sec)');
    pp_thresh = 10^(30/20);
    add_analysis_callbacks;
    rc = evalin('base','ref_channel;');
    if(rc == 1)
        pp = single(abs(v1) > pp_thresh);
    else
        pp = single(abs(v2) > pp_thresh);
    end
    pp(pp == 0) = nan;
    
    figure;
    subplot(2,1,1);
    if(rc == 1)
        plot(ts*(0:length(v1)-1),pp.*[angle(v1)],'.-');
    else
        plot(ts*(0:length(v1)-1),[angle(v1)],'.-');
    end
    grid on;
    %xlabel('Frequency (MHz)');
    %ylabel('Amplitude (20log_{10})');
    title('Phase (rad) [CH1]');ah = [ah gca];
    xlabel('Time (sec)');

    subplot(2,1,2);
    if(rc == 1)
        plot(ts*(0:length(v1)-1),[angle(v2)],'.-');
    else
        plot(ts*(0:length(v2)-1),pp.*[angle(v2)],'.-');
    end
        grid on;
    %xlabel('Frequency (MHz)');
    %ylabel('Amplitude (20log_{10})');
    title('Phase (rad) [CH2]');ah = [ah gca];
    
    xlabel('Time (sec)');
    
    add_analysis_callbacks;
    
    figure;
    plot(ts*(0:length(v1)-1),[pp.*angle(v1./v2)]./(2*pi),'.-');
    grid on;
    xlabel('Time (sec)');
    ylabel('Cycles');
    title('Rel Phase (rad) [CH2]');ah = [ah gca];
    linkaxes(ah,'x');
    xlabel('Time (sec)');
    
    add_analysis_callbacks;
    
%end


end

