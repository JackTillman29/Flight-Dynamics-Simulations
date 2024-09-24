close all;
clear;
clc;

addpath('z:\sawmillkd\MATLAB\DigitalFilters');
addpath('z:\sawmillkd\MATLAB\DSP');
addpath('z:\sawmillkd\MATLAB\Windows');


% Load file
BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2204 D20 Mod 1 NBG V4_1_Alt_T STIM.dat';
%BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2206 D20 WBG17 Baseline V4_1_Alt_T STIM.dat';





[fp_path,fp_file,fp_ext]=fileparts(BSAG_FILE);
BSAG_HEADER = [fp_path '\' fp_file '.txt'];
headerData = ReadBSAGHeader(BSAG_HEADER);
dt_data  = 1/(headerData.DATA_SAMPLE_RATE_MHZ*1e6);

flags.showFirstTrigger = 0;

previewHop = 1; % number of data frames to hop (5 samples per frame)
nPacketsPerRead = 2000000;
nReadsInFile    =(headerData.SIGNAL_NUM_SAMPLES/5)/nPacketsPerRead/previewHop;

thresholds.pulserise = 18; % counts magnitude
thresholds.pulsefall = 8; % counts magnitude
thresholds.hardzero = 160; % counts magnitude for pulses
triggering.presamples  = 1000/previewHop;
triggering.postsamples = (5.25e-3/(previewHop*5))/dt_data;

% pulse level triggering
thresholds.pulserise_pulse = 200; % counts magnitude


fir.burst=firFilter(0,20000,1/dt_data,256*16,'blackman(x)');

%=========================================================================

[fp_path,fp_file,fp_ext]=fileparts(BSAG_FILE);
BSAG_HEADER = [fp_path '\' fp_file '.txt'];
headerData = ReadBSAGHeader(BSAG_HEADER);
fid = fopen(BSAG_FILE,'rb');

fseek(fid,0,'eof');
nBytes = ftell(fid);

if(nBytes ~= headerData.SIGNAL_NUM_BYTES)
    error('Header bytes does not match file size');
end




% Work on loop to read and parse files
% output file
BSAG_BURSTFILE = [fp_file '_burst.txt'];
fid_burst = fopen(BSAG_BURSTFILE,'w');

dt_preview = dt_data * 5 * previewHop;



startingPacket = 1;


firstPacket = startingPacket;
tic;
for k = 1 : nReadsInFile
    lastPacket = firstPacket - 1 + nPacketsPerRead*previewHop;
    [data,time,packet_num]=ReadBSAG12x5(fid,headerData,'preview',previewHop,firstPacket,nPacketsPerRead);
    triggers = TriggerPoints(fir.burst,[thresholds.pulserise thresholds.pulsefall], ...
        triggering.presamples,triggering.postsamples,single(data.m'),flags.showFirstTrigger);
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
        disp(['             [data,time,packet_num]=ReadBSAG12x5(fid,headerData,''full'',[],' num2str(triggers.lead_packet(itrig)+firstPacket-1) ',' num2str(triggers.trail_packet(itrig)-triggers.lead_packet(itrig)) ');']);
        if(flags.showFirstTrigger)
            return;
            %%
        end
    end
    disp(['Rate: ' num2str(lastPacket/toc) ' pps']);
    firstPacket = lastPacket + 1;
    
    
end
fclose all;
%%
sumdat = load(BSAG_BURSTFILE);
figure;
plot(sumdat(:,3),'.');
ylabel('Pulsewidth(ns)');
title('Pulsewidth Summary');

figure;
plot(sumdat(:,4),'.');
ylabel('PRF (Hz)');
title('PRF Summary');

sumdat(3240,:)