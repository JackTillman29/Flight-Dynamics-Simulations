close all;
fclose all;
clear;
clc;

addpath('z:\sawmillkd\MATLAB\DigitalFilters');
addpath('z:\sawmillkd\MATLAB\DSP');

% Load file
BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2193 D20 NBG Baseline V4_1 STIM.dat';
BSAG_FILE_TG  = 'E:\Data from 20150622\BSAG TG\T2193 D20 NBG Baseline V4_1 TG.dat';

[fp_path,fp_file,fp_ext]=fileparts(BSAG_FILE);
[fp_path_TG,fp_file_TG,fp_ext_TG]=fileparts(BSAG_FILE_TG);

BSAG_HEADER = [fp_path '\' fp_file '.txt'];
BSAG_HEADER_TG = [fp_path_TG '\' fp_file_TG '.txt'];

headerData = ReadBSAGHeader(BSAG_HEADER);
headerData_TG = ReadBSAGHeader(BSAG_HEADER_TG);

% output file



%%
fid = fopen(BSAG_FILE,'rb');
fid_TG = fopen(BSAG_FILE_TG,'rb');

dt_data  = 1/(headerData.DATA_SAMPLE_RATE_MHZ*1e6);
maxEnvelope = sqrt(headerData.LOADING_MAX_VALUE.^2 + headerData.LOADING_MIN_VALUE.^2);
maxEnvelope_TG = sqrt(headerData_TG.LOADING_MAX_VALUE.^2 + headerData_TG.LOADING_MIN_VALUE.^2);
skipval = 0;



% firstSample = 1053871; %422957413; 
% lastSample  = 1124846; %426759794;

%firstSample = 469014524; % 83.5khz
%lastSample  = 469085499; % 

%firstSample = 424786536 ; % 81.9khz
%lastSample  = 424857511 ; % 

firstSample = 423359818 ; % 95.8khz
lastSample  = 423430793 ; % 
 

nSamples = lastSample - firstSample + 1;

fseek(fid,(firstSample-1)*headerData.DATA_BYTES_PER_UNIT,'bof');
fseek(fid_TG,(firstSample-1)*headerData_TG.DATA_BYTES_PER_UNIT,'bof');

% for burst
thresholds.hardzero = 100000;
thresholds.pulserise = 6e6;
thresholds.maxysize  = 216477760;
triggering.presamples  = 100;
triggering.postsamples = 5.4e-3/dt_data;
% for pulse
thresholds.pulserise_indiv = 2e7;
triggering.pulserise_indiv_pre   = 100;
triggering.pulserise_indiv_post  = floor(10e-6 / dt_data);



fir.burst=firFilter(0,20000,1/dt_data,256*16,'blackman(x)');
%fir.pulse=firFilter(0,1e6,1/dt_data,256*16,'blackman(x)');applyFilter(fir.pulse,single(xm2')); not needed

figure(1);
h.burstplot_data = plot(0,0);
h.burstplot_axes = gca;
set(h.burstplot_axes,'YLim',[0 thresholds.maxysize]);
tic;
%for k = 1 : 1
    lastSample = firstSample-1+nSamples;
    %disp([k firstSample lastSample 100*firstSample / headerData.SIGNAL_NUM_SAMPLES]);
    % read counts as signed int16
    x   = fread(fid, 2*nSamples, '*int16', skipval, 'ieee-be');
    x_TG = fread(fid_TG, 2*nSamples, '*int16', skipval, 'ieee-be');
    % compute iq magnitude squared
    xm2 = int32(x(1:2:end)).^2+int32(x(2:2:end)).^2;
    xm2_TG = int32(x_TG(1:2:end)).^2+int32(x_TG(2:2:end)).^2;

    % Triggers for bursts
    triggers = TriggerPoints(fir.burst,thresholds.pulserise,triggering.presamples,triggering.postsamples,single(xm2'),1);
    triggers.mid = triggers.lead + floor(0.5*(triggers.trail - triggers.lead)); 
    
    % Read pulsewidth & prf data for this burst
    if(~isempty(triggers.lead))
    for itrig = 1 : length(triggers.lead)
        [pw,prf,prf2,pwr] = GetWaveformProps(thresholds.pulserise,thresholds.hardzero,xm2(triggers.lead(itrig):triggers.trail(itrig)),dt_data);
        ydata = xm2(triggers.lead(itrig):triggers.trail(itrig));
        set(h.burstplot_data,'XData',(0:(length(ydata)-1)),'YData',ydata);
        %drawnow;
        %fprintf(fid_burst,'%d %d %d %d %d %d\n',[triggers.lead(itrig)+firstSample-1 triggers.trail(itrig)+firstSample-1 floor(1e9*pw) floor(prf) floor(prf2) floor(pwr)]);
        fprintf(         'PCT=%3d Start: %12d %8dns %6dHz %6dHz Amp: %6d\n',[floor(100*firstSample/headerData.SIGNAL_NUM_SAMPLES) firstSample floor(1e9*pw) floor(prf) floor(prf2) floor(pwr)]);
        h.text(itrig) = text(0,0,'test');
        set(h.text(itrig),'Position',[triggers.mid(itrig) 0],'Rotation',90,'Color','r','FontWeight','bold');
        set(h.text(itrig),'String',num2str([prf2/1000],'%8.3f kHz'));
        
    end
    else
        disp('here');
        triggers_pulse = TriggerPoints(xm2',thresholds.pulserise_indiv,triggering.pulserise_indiv_pre,triggering.pulserise_indiv_post,single(xm2'),1);
        triggers_pulse.mid = triggers_pulse.lead + floor(0.5*(triggers_pulse.trail - triggers_pulse.lead));
        imData = zeros(triggers_pulse.trail-triggers_pulse.lead+1,length(triggers_pulse.lead),'single')';
        imData_TG = 0*imData;
        
        % generate timing jitter in data post-test
        
        
        for imRow = 1 : size(imData,1)
            idx_start = 2*triggers_pulse.lead(imRow)-1;
            idx_stop  = 2*triggers_pulse.trail(imRow)-1;
            idx_offset = GenerateTimingJitter(size(imData,1),floor(.5e-6/dt_data));
            
            imData(imRow,:) = abs(single(x(idx_start:2:idx_stop)) + 1i * single(x((idx_start+1):2:(idx_stop+1))));
            imData_TG(imRow,:) = abs(single(x_TG((idx_start+idx_offset):2:(idx_stop+idx_offset))) + ...
                                1i * single(x_TG((idx_start+1+idx_offset):2:(idx_stop+1+idx_offset))));
        end
        
        xVec = (triggers_pulse.lead(1):triggers_pulse.trail(1))-triggers_pulse.mid(1) + (triggering.pulserise_indiv_post-triggering.pulserise_indiv_pre)/4;
        yVec = triggers_pulse.mid;
        xVec = xVec * dt_data * 1e6;
        yVec = yVec * dt_data * 1e6;
        %%
        figure;
        subplot(2,2,1);
        imagesc(xVec,yVec,imData);
        title('Stimulus');
        grid on;
        xlabel('PRI Time (\mus)');ylabel('Burst Time (\mus)');colorbar;ah = gca;
        
        subplot(2,2,2);
        imagesc(xVec,yVec,imData_TG);
        title('TG Response');
        grid on;
        xlabel('PRI Time (\mus)');ylabel('Burst Time (\mus)');colorbar;ah = [ah gca];
        linkaxes(ah,'xy');
        
        subplot(2,2,[3 4]);
        plot(xVec,sum(imData,1)./size(imData,1),xVec,sum(imData_TG,1)./size(imData_TG,1));
        grid on;
        xlabel({'PRI Time (\mus)',fp_file});ylabel('A/D Counts');
        title(['Average waveform over ' num2str(size(imData_TG,1)) ' samples.']);
        
    end
    
    firstSample=lastSample+1;
%end
elapsedTime = toc;
sampPerSec = firstSample / elapsedTime;
disp(['Processing Rate: ' num2str(sampPerSec) ' samples / sec.']);
fclose(fid);
fclose(fid_TG);
%fclose(fid_burst);