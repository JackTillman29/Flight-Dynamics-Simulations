close all;
fclose all;
clear;
clc;

addpath('z:\sawmillkd\MATLAB\DigitalFilters');
addpath('z:\sawmillkd\MATLAB\DSP');
addpath('z:\sawmillkd\MATLAB\Windows');
addpath('z:\sawmillkd\MATLAB\Printing');

% Load file
%BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2209 S20 NBG Baseline V4_1_Phase Correction STIM.dat';
%BSAG_FILE_TG  = 'E:\Data from 20150622\BSAG TG\T2209 S20 NB Baseline V4_1_Phase Correction TG.dat';
BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2204 D20 Mod 1 NBG V4_1_Alt_T STIM.dat';
BSAG_FILE_TG  = 'E:\Data from 20150622\BSAG TG\T2204 D20 Mod1 NBG V4_1_Alt_T TG.dat';


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

                 
firstPacket = 1431316256;
lastPacket  = 1431537756 ;


    
 

nPackets = lastPacket - firstPacket + 1;

fseek(fid,(firstPacket-1)*headerData.DATA_BYTES_PER_UNIT,'bof');
fseek(fid_TG,(firstPacket-1)*headerData_TG.DATA_BYTES_PER_UNIT,'bof');

% for pulse
thresholds.pulserise_indiv = 75;
triggering.pulserise_indiv_pre   = floor(10e-6 / dt_data);
triggering.pulserise_indiv_post  = floor(10e-6 / dt_data);

fir.burst=firFilter(0,80000,1/dt_data,32,'blackman(x)');
%fir.pulse=firFilter(0,1e6,1/dt_data,256*16,'blackman(x)');applyFilter(fir.pulse,single(xm2')); not needed

figure(1);
h.burstplot_data = plot(0,0);
h.burstplot_axes = gca;
%set(h.burstplot_axes,'YLim',[0 thresholds.maxysize]);
tic;
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
    %TG
    [data,time,packet_num]=ReadBSAG12x5(fid_TG,headerData,'full',[],firstPacket,nPackets);
    xm2_TG = data.m;

    % Triggers for bursts
    triggers_pulse = TriggerPoints(xm2',thresholds.pulserise_indiv,triggering.pulserise_indiv_pre,triggering.pulserise_indiv_post,single(xm2),1);
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
        xlabel({'PRI Time (\mus)',strrep(fp_file,'_','\_')});ylabel('A/D Counts');
        title(['Average waveform over ' num2str(size(imData_TG,1)) ' samples.']);
        
    

%end
%elapsedTime = toc;
%sampPerSec = firstSample / elapsedTime;
%disp(['Processing Rate: ' num2str(sampPerSec) ' samples / sec.']);
fclose(fid);
fclose(fid_TG);
%fclose(fid_burst);
%%
pp = PrepForPrint();
PrepForPrint(gcf,pp);
grid on;
set(gca,'XColor',[0 0 0]+0.7,'YColor',[0 0 0]+0.7);
fo=findobj('Parent',gca);
set(fo(2),'Color',[.9 .9 1]);
kprint(gcf,'bsag_dyn20_mod1_95p67khz.png',180);