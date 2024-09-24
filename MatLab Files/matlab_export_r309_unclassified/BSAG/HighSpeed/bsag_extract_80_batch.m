close all;
fclose all;
clear;
clc;

addpath('z:\sawmillkd\MATLAB\DigitalFilters');
addpath('z:\sawmillkd\MATLAB\DSP');
addpath('z:\sawmillkd\MATLAB\Windows');
addpath('z:\sawmillkd\MATLAB\Printing');
pp = PrepForPrint();
% Load file
%BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2209 S20 NBG Baseline V4_1_Phase Correction STIM.dat';
%BSAG_FILE_TG  = 'E:\Data from 20150622\BSAG TG\T2209 S20 NB Baseline V4_1_Phase Correction TG.dat';

%---------------------------------
BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2204 D20 Mod 1 NBG V4_1_Alt_T STIM.dat';
BSAG_FILE_TG  = 'E:\Data from 20150622\BSAG TG\T2204 D20 Mod1 NBG V4_1_Alt_T TG.dat';
prf_batch = [94180 82000 83713 95673];                 
firstPacket_arr = [98199031 1504364231 1426751106 1431316256];
lastPacket_arr  = [98420531 1504585731 1426972606 1431537756];

sumdat = load('T2204 D20 Mod 1 NBG V4_1_Alt_T STIM_burst.txt');
sumdat(:,4) = floor(sumdat(:,4)/100)*100;
prf_batch = sumdat(:,4);
firstPacket_arr = sumdat(:,1);
lastPacket_arr = sumdat(:,2);
fidbw = fopen('avgdata.bin','wb');

%-----------------------------------
%BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2206 D20 WBG17 Baseline V4_1_Alt_T STIM.dat';
%BSAG_FILE_TG  = 'E:\Data from 20150622\BSAG TG\T2206 D20 WBG17 Baseline V4_1_Alt_T TG.dat';
%prf_batch = [82000 94179 95670 83710];                 
%firstPacket_arr = [1319446400 1319680458 1328577274 1337708273];
%lastPacket_arr  = [1319667900 1319901958 1328798774 1337929773];
%-----------------------------------


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
%maxEnvelope = sqrt(headerData.LOADING_MAX_VALUE.^2 + headerData.LOADING_MIN_VALUE.^2);
%maxEnvelope_TG = sqrt(headerData_TG.LOADING_MAX_VALUE.^2 + headerData_TG.LOADING_MIN_VALUE.^2);
skipval = 0;



% firstSample = 1053871; %422957413; 
% lastSample  = 1124846; %426759794;

%firstSample = 469014524; % 83.5khz
%lastSample  = 469085499; % 

%firstSample = 424786536 ; % 81.9khz
%lastSample  = 424857511 ; % 

%firstPacket = 1654951 ;
%lastPacket  = 1876451 ; % 'E:\Data from 20150622\BSAG Stim\T2209 S20 NBG Baseline V4_1_Phase Correction STIM.dat';


pri_batch = 1./prf_batch;

    
    
for kbatch = 4020 : length(prf_batch) 
firstPacket = firstPacket_arr(kbatch);
lastPacket  = lastPacket_arr(kbatch);

nPackets = lastPacket - firstPacket + 1;

fseek(fid,(firstPacket-1)*headerData.DATA_BYTES_PER_UNIT,'bof');
fseek(fid_TG,(firstPacket-1)*headerData_TG.DATA_BYTES_PER_UNIT,'bof');

% for pulse
thresholds.pulserise_indiv = 75;
triggering.pulserise_indiv_pre   = 10; %floor(10e-6 / dt_data);
triggering.pulserise_indiv_post  = floor(pri_batch(kbatch) / dt_data);

fir.burst=firFilter(0,80000,1/dt_data,32,'blackman(x)');
%fir.pulse=firFilter(0,1e6,1/dt_data,256*16,'blackman(x)');applyFilter(fir.pulse,single(xm2')); not needed

%figure(1);
%h.burstplot_data = plot(0,0);
%h.burstplot_axes = gca;
%set(h.burstplot_axes,'YLim',[0 thresholds.maxysize]);
%tic;
%for k = 1 : 1
  %  lastPacket = firstPacket-1+nPackets
    %disp([k firstSample lastSample 100*firstSample / headerData.SIGNAL_NUM_SAMPLES]);
    % read counts as signed int16
   % x   = fread(fid, 2*nSamples, '*int16', skipval, 'ieee-be');
    %x_TG = fread(fid_TG, 2*nSamples, '*int16', skipval, 'ieee-be');
    % compute iq magnitude squared
   % xm2 = int32(x(1:2:end)).^2+int32(x(2:2:end)).^2;
   % xm2_TG = int32(x_TG(1:2:end)).^2+int32(x_TG(2:2:end)).^2;
    
    %STIM
    [data,time,packet_num]=ReadBSAG12x5(fid,headerData,'full',[],firstPacket,nPackets);
    xm2 = data.m;
    close(gcf);
    %TG
    [data,time,packet_num]=ReadBSAG12x5(fid_TG,headerData,'full',[],firstPacket,nPackets);
    xm2_TG = data.m;
    close(gcf);

    % Triggers for bursts
    triggers_pulse = TriggerPoints(xm2',thresholds.pulserise_indiv,triggering.pulserise_indiv_pre,triggering.pulserise_indiv_post,single(xm2),0);
    triggers_pulse.mid = triggers_pulse.lead + floor(0.5*(triggers_pulse.trail - triggers_pulse.lead));
    imData = zeros(triggers_pulse.trail-triggers_pulse.lead+1,length(triggers_pulse.lead),'single')';
    imData_TG = 0*imData;
        
        % generate timing jitter in data post-test
        
         for imRow = 1 : size(imData,1)
            idx_start = triggers_pulse.lead(imRow);
            idx_stop  = triggers_pulse.trail(imRow);
            idx_offset = 0;
            %idx_offset = GenerateTimingJitter(size(imData,1),floor(0.5e-6/dt_data));
            
            imData(imRow,:) = single(xm2(idx_start:idx_stop));
            imData_TG(imRow,:) = single(xm2_TG((idx_start+idx_offset):(idx_stop+idx_offset)));
        end
        
        xVec = (triggers_pulse.lead(1):triggers_pulse.trail(1))-triggers_pulse.mid(1) + (triggering.pulserise_indiv_post-triggering.pulserise_indiv_pre)/4;
        xVec = xVec - xVec(1);
        yVec = triggers_pulse.mid;
        xVec = xVec * dt_data * 1e6;
        yVec = yVec * dt_data * 1e6;
        %%
        figure;
        set(gcf,'Name',[num2str(prf_batch(kbatch)) 'Hz Image']);
        subplot(1,2,1);
        imagesc(xVec,yVec,imData);
        title('Stimulus');
        grid on;
        xlabel('PRI Time (\mus)');ylabel('Burst Time (\mus)');colorbar;ah = gca;
        
        subplot(1,2,2);
        imagesc(xVec,yVec,imData_TG);
        title('TG Response');
        grid on;
        xlabel('PRI Time (\mus)');ylabel('Burst Time (\mus)');colorbar;ah = [ah gca];
        linkaxes(ah,'xy');
        set(gcf,'Position',[1304 508 1862 420]);
        PrepForPrint(gcf,pp);
        kprint(gcf,[fp_file ' ' num2str(prf_batch(kbatch)) 'Hz Image ' num2str(firstPacket) ' ' num2str(lastPacket) '.png'],96);
        
        figure;
        set(gcf,'Name',[num2str(prf_batch(kbatch)) 'Hz Average']);
        plot(xVec,sum(imData,1)./size(imData,1),xVec,sum(imData_TG,1)./size(imData_TG,1));
        grid on;
        xlabel({'PRI Time (\mus)',strrep(fp_file,'_','\_')});ylabel('A/D Counts');
        title(['Average waveform over ' num2str(size(imData_TG,1)) ' PRIs']);
        set(gcf,'Position',[1304 508 1862 420]);
        
        fwrite(fidbw,[size(imData,2) sum(imData,1)./size(imData,1) sum(imData_TG,1)./size(imData_TG,1)],'single');

        PrepForPrint(gcf,pp);
        kprint(gcf,[fp_file ' ' num2str(prf_batch(kbatch)) 'Hz Average ' num2str(firstPacket) ' ' num2str(lastPacket) '.png'],96);
          
        % grid on;
        %set(gca,'XColor',[0 0 0]+0.7,'YColor',[0 0 0]+0.7);
        %fo=findobj('Parent',gca);
        %set(fo(2),'Color',[.9 .9 1]);
        
        x = single(xm2);
        y = single(xm2_TG);
        figure;
        xf=applyFilter(fir.burst,x(:)')>50;
        yf=applyFilter(fir.burst,y(:)')>50;
        %plot(time,[x(:) y(:)]);
        dc = 10;
        plot(1e6*time(1:dc:end),[xf(1:dc:end)' 1.1+yf(1:dc:end)']);
        title(['Technique Timing Diagram: ' num2str(prf_batch(kbatch)) 'Hz']);
        set(gca,'YTick',[0.5 1.6],'YTickLabel',{'Threat','TG'});
        set(gcf,'Name',[num2str(prf_batch(kbatch)) 'Hz Rolling']);
        set(gcf,'Position',[1304 508 1862 420]);
        xlabel('Time (\mus)');
        PrepForPrint(gcf,pp);
        kprint(gcf,[fp_file ' ' num2str(prf_batch(kbatch)) 'Hz Timing Full ' num2str(firstPacket) ' ' num2str(lastPacket) '.png'],96);
        xlim([300 800]);
        kprint(gcf,[fp_file ' ' num2str(prf_batch(kbatch)) 'Hz Timing Zoom ' num2str(firstPacket) ' ' num2str(lastPacket) '.png'],96);
        close all;
end
%end
%elapsedTime = toc;
%sampPerSec = firstSample / elapsedTime;
%disp(['Processing Rate: ' num2str(sampPerSec) ' samples / sec.']);
fclose(fid);
fclose(fid_TG);
fclose(fidbw);
%fclose(fid_burst);


%kprint(gcf,'bsag_dyn20_mod1_95p67khz.png',180);