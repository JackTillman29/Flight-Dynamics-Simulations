function DS1054Z_RealTimeCh1(obj)
    %% This function requires that the SCPI I/F be "closed"
    fclose(obj);
    obj.InputBufferSize = 16384;
    warning('off','instrument:fread:unsuccessfulRead');
    ichan = 1;
    
    fopen(obj);
    fprintf(obj,':WAV:MODE MAX')
    fprintf(obj,':WAV:FORM BYTE')
    fprintf(obj,[':WAV:SOUR CHAN' num2str(ichan)])
    

    
    init = 1;
    fps_counter = 1;
    
    while(1)
        if(fps_counter == 1)
            tic;
            fprintf(obj,':WAV:XINC?');
            xinc = str2num(char(fread(obj)'));
            fprintf(obj,':WAV:XREF?');
            xref = str2num(char(fread(obj)'));
            fprintf(obj,':WAV:XOR?');
            xorig = str2num(char(fread(obj)'));
            Fs = 1 / xinc;

            fprintf(obj,':WAV:YINC?');
            yinc = str2num(char(fread(obj)'));
            fprintf(obj,':WAV:YREF?');
            yref = str2num(char(fread(obj)'));
            fprintf(obj,':WAV:YOR?');
            yorig = str2num(char(fread(obj)'));
        end
        fps_counter = fps_counter + 1;
        fprintf(obj,':WAV:DATA?');
        [data,len]=fread(obj);
        V = (single(data(12:end-1))-yref)*yinc;
        fV = 20*log10(abs(fft(hamming(length(V)).*V)));
        fV = fV(1:floor(length(fV)/2));
        fVx = (0:length(fV)-1)./(length(fV));
        fVx = fVx ./ xinc ./ 2;
        
        if(init)
            hf = figure;
            %hL = plot(V);
            hL = plot(fV);
            hT = title(['Fs = ' num2str(Fs)]);
            init = 0;
        else
            %set(hL,'YData',V);
            set(hL,'YData',fV);
            set(hT,'String',['Fs = ' num2str(Fs)]);
            drawnow;
        end
        %disp(len);
        if(fps_counter == 10)
            frame_time = toc;
            frame_time = 0.1*frame_time;
            fps_counter = 1;
            disp(['fps = ' num2str(1/frame_time)]);
            set(hL,'XData',fVx);
        end
    end
    
end