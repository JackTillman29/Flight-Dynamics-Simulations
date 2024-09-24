close all;
clear;
clc;

addpath('z:\sawmillkd\MATLAB\DigitalFilters');
addpath('z:\sawmillkd\MATLAB\DSP');
addpath('z:\sawmillkd\MATLAB\Windows');


% Load file
%BSAG_FILE = 'E:\Data from 20150622\BSAG Stim\T2209 S20 NBG Baseline V4_1_Phase Correction STIM.dat';
%BSAG_FILE_TG = 'E:\Data from 20150622\BSAG TG\T2209 S20 NB Baseline V4_1_Phase Correction TG.dat';

%BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2218 D20 NBG Baseline V4_1_Mod_4_LG STIM.dat';
BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2204 D20 Mod 1 NBG V4_1_Alt_T STIM.dat';

[fp_path,fp_file,fp_ext]=fileparts(BSAG_FILE);
BSAG_HEADER = [fp_path '\' fp_file '.txt'];
avifilename = [fp_path '\' fp_file '.avi'];
headerData = ReadBSAGHeader(BSAG_HEADER);
fid = fopen(BSAG_FILE,'rb');

fseek(fid,0,'eof');
nBytes = ftell(fid);

if(nBytes ~= headerData.SIGNAL_NUM_BYTES)
    error('Header bytes does not match file size');
end
dt_data  = 1/(headerData.DATA_SAMPLE_RATE_MHZ*1e6);


if(0)
%% Checks out so far... read in some real data
startingPacket = 1;
nPackets = 1000000;

previewMode = 1;
if(~previewMode)
    [data,time,sample_num]=ReadBSAG12x5(fid,headerData,'full'   ,[],startingPacket,nPackets);
else
    [data,time,packet_num]=ReadBSAG12x5(fid,headerData,'preview',1,startingPacket,nPackets);
end
return;
end

if(0)
    %% Testing a full read w/ image processing
    fid = fopen(BSAG_FILE,'rb');
    fid_TG = fopen(BSAG_FILE_TG,'rb');
    startingPacket = 1654951    ;
    endingPacket   = 1876451  ;
    thresholds.pulserise_indiv = 100;
    
    nPackets       = endingPacket - startingPacket + 1;
    [data,time,sample_num]=ReadBSAG12x5(fid,headerData,'full',[],startingPacket,nPackets);
    ah = gca;
    [data_TG]=ReadBSAG12x5(fid_TG,headerData,'full',[],startingPacket,nPackets);
    ah = [ah gca];
    linkaxes(ah,'xy');
    triggering.pulserise_indiv_pre =1100;
    triggering.pulserise_indiv_post=1100;
    triggers_pulse = TriggerPoints(data.m,[thresholds.pulserise_indiv 50],triggering.pulserise_indiv_pre,triggering.pulserise_indiv_post,data.m,1);
    triggers_pulse.mid = triggers_pulse.lead + floor(0.5*(triggers_pulse.trail - triggers_pulse.lead));
    imData = zeros(triggers_pulse.trail-triggers_pulse.lead+1,length(triggers_pulse.lead),'single')';
    imData_TG = 0*imData;
    
    for imRow = 1 : size(imData,1)
            idx_start = triggers_pulse.lead(imRow)-1;
            idx_stop  = triggers_pulse.trail(imRow)-1;
            %idx_offset = GenerateTimingJitter(size(imData,1),floor(.5e-6/dt_data));
            
            imData(imRow,:) = single(data.m(idx_start:idx_stop));
            %imData_TG(imRow,:) = abs(single(x_TG((idx_start+idx_offset):2:(idx_stop+idx_offset))) + ...
            %                    1i * single(x_TG((idx_start+1+idx_offset):2:(idx_stop+1+idx_offset))));
            imData_TG(imRow,:) = single(data_TG.m(idx_start:idx_stop));
    end

    xVec = (triggers_pulse.lead(1):triggers_pulse.trail(1))-triggers_pulse.mid(1) + (triggering.pulserise_indiv_post-triggering.pulserise_indiv_pre)/4;
    yVec = triggers_pulse.mid;
    xVec = xVec * dt_data * 1e6;
    yVec = yVec * dt_data * 1e6;
    
    
    
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

previewHop = 2; % number of data frames to hop (5 samples per frame)

% Work on loop to read and parse files
% output file
BSAG_BURSTFILE = [fp_path '\' fp_file '_burst.txt'];
fid_burst = fopen(BSAG_BURSTFILE,'w');

dt_preview = dt_data * 5 * previewHop;

nPacketsPerRead = 1000000;
nReadsInFile    =(headerData.SIGNAL_NUM_SAMPLES/5)/nPacketsPerRead/previewHop;

thresholds.pulserise = 18; % counts magnitude
thresholds.pulsefall = 8; % counts magnitude
thresholds.hardzero = 160; % counts magnitude for pulses
triggering.presamples  = 1000/previewHop;
triggering.postsamples = (5.25e-3/(previewHop*5))/dt_data;

% pulse level triggering
thresholds.pulserise_pulse = 200; % counts magnitude


fir.burst=firFilter(0,20000,1/dt_data,256*16,'blackman(x)');


startingPacket = 1;


firstPacket = startingPacket;
tic;
for k = 1 : nReadsInFile
    lastPacket = firstPacket - 1 + nPacketsPerRead*previewHop;
    [data,time,packet_num]=ReadBSAG12x5(fid,headerData,'preview',previewHop,firstPacket,nPacketsPerRead);
    triggers = TriggerPoints(fir.burst,[thresholds.pulserise thresholds.pulsefall], ...
        triggering.presamples,triggering.postsamples,single(data.m'),0);
    [triggers.lead_packet,triggers.lead_sample] = P2F(triggers.lead,previewHop);
    [triggers.trail_packet,triggers.trail_sample] = P2F(triggers.trail,previewHop);
    
    % Read pulsewidth & prf data for this burst
    for itrig = 1 : length(triggers.lead)
        ydata = data.m(triggers.lead(itrig):triggers.trail(itrig));
        [pw,prf,prf2,pwr] = GetWaveformProps(thresholds.pulserise_pulse,thresholds.hardzero,ydata,dt_preview);
        fprintf(fid_burst,'%d %d %d %d %d %d\n',[ ...
            triggers.lead_packet(itrig)+firstPacket-1 ...
            triggers.trail_packet(itrig)+firstPacket-1 ...
            floor(1e9*pw) floor(prf) floor(prf2) floor(pwr)]);
        fprintf(         'PCT=%3d Start: %12d %8dns %6dHz %6dHz Amp: %6d\n',[floor(100*firstPacket/(headerData.SIGNAL_NUM_SAMPLES/5)) firstPacket floor(1e9*pw) floor(prf) floor(prf2) floor(pwr)]);
%         figure;
%         plot(ydata);
%         pause;
%         return;
    end
    disp(['Rate: ' num2str(lastPacket/toc) ' pps']);
    firstPacket = lastPacket + 1;
    
    
end
fclose all;