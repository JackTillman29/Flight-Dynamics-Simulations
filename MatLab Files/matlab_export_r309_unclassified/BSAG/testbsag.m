close all;
clear all;
clc;

% Load file
BSAG_FILE  = 'E:\Data from 20150622\BSAG Stim\T2188 S20 NBG Baseline V4.1 Stim.dat';
BSAG_FILE2 = 'E:\Data from 20150622\BSAG TG\T2188 S20 NBG Baseline V4.1 TG.dat';

[fp_path,fp_file,fp_ext]=fileparts(BSAG_FILE);
BSAG_HEADER = [fp_path '\' fp_file '.txt'];
avifilename = [fp_path '\' fp_file '.avi'];
headerData = ReadBSAGHeader(BSAG_HEADER);
fid = fopen(BSAG_FILE,'rb');

[fp_path2,fp_file2,fp_ext2]=fileparts(BSAG_FILE2);
BSAG_HEADER2 = [fp_path2 '\' fp_file2 '.txt'];
headerData2 = ReadBSAGHeader(BSAG_HEADER2);
fid2 = fopen(BSAG_FILE2,'rb');

fseek(fid2,0,'eof');
nBytes2 = ftell(fid2);
frewind(fid2);
if(nBytes2 ~= headerData2.SIGNAL_NUM_BYTES)
    error('Header bytes does not match file size');
end




%% Purpose of this file is to generate avi files for BSAG data to see what 
% is in the data
vidrec = 0;
fps = 10;
dt_frame = 1 / fps;
dt_data  = 1/(headerData.DATA_SAMPLE_RATE_MHZ*1e6);
nSamples = floor(dt_frame / dt_data);
nFrames = floor(headerData.SIGNAL_NUM_SAMPLES / nSamples);

maxEnvelope = sqrt(headerData.LOADING_MAX_VALUE.^2 + headerData.LOADING_MIN_VALUE.^2);

skipval = 0;
figure;
set(1,'Position',[13         282        1210         636]);
set(1,'PaperPositionMode','auto');

firstSample = 1;



for k = 1 : 2 %nFrames
    lastSample = firstSample-1+nSamples;
    disp([k firstSample lastSample]);
    x = fread(fid, 2*nSamples, '*int16', skipval, 'ieee-be');
    xq = single(x(1:2:end))+1i*single(x(2:2:end));
    x2 = fread(fid2, 2*nSamples, '*int16', skipval, 'ieee-be');
    xq2 = single(x2(1:2:end))+1i*single(x2(2:2:end));
    
    
    if(k==1)
        figure(1);
        subplot(2,1,1);
        hData = plot(firstSample:lastSample,abs(xq));
        hAxis = gca;
        ht=title(dt_data*firstSample);
        subplot(2,1,2);
        hData2 = plot(firstSample:lastSample,abs(xq2));
        hAxis2 = gca;
        if(vidrec)
            vidObj = VideoWriter(avifilename);
            vidObj.FrameRate = 30.0;
            vidObj.Quality = 95;
            open(vidObj);
        end
        set([hAxis hAxis2],'YLim',[0 maxEnvelope]);
        
        this_xtick = get(hAxis,'XTick');
        set([hAxis hAxis2],'XTickLabel',num2str(this_xtick'*dt_data*1000,'%8.2f'));
        set(hAxis2,'Position',[0.025 0+0.025 1-2*0.025 .5-2*0.025]);
        set(hAxis,'Position',[0.025 .5 1-2*0.025 .5-2*0.025]);
        set(gcf,'Color','k');
        set([hAxis hAxis2],'XColor','w','YColor','w','Color','k','XGrid','on','YGrid','on');
        set([hData hData2],'Color','g');
        drawnow;
    else
        set(hData,'YData',abs(xq),'XData',firstSample:lastSample);
        set(hData2,'YData',abs(xq2),'XData',firstSample:lastSample);
        
        set(ht,'String',dt_data*firstSample);
        this_xtick = get(hAxis,'XTick');
        set([hAxis hAxis2],'XTickLabel',num2str(this_xtick'*dt_data*1000,'%8.2f'));
        drawnow;
        linkaxes([hAxis hAxis2],'x');
    end
  
    %print('-dpng','-f1','-r180',['Frame_' num2str(k) '.png']);
    if(vidrec)
        writeVideo(vidObj,hardcopy(gcf,'-Dopengl'));
    end
    

    firstSample=lastSample+1;
end
fclose(fid);
if(vidrec)
    close(vidObj);
end

%% Filter 
%u = xq2(1:4e5);
close all;
addpath('z:\sawmillkd\MATLAB\PulsesAndWaveforms');
addpath('z:\sawmillkd\MATLAB\Utilities');
u = xq(75374:148899);
usig = xq2(75374:148899);


dt = 1/(headerData.DATA_SAMPLE_RATE_MHZ*1e6);
tu = (0:(length(u)-1))'*dt;

fif_1 = 12.055e6;
%f2ndlo = 48.220e6;

dt2 = 1/(10*fif_1);
tu2 = single(0:dt2:tu(end));
u2 = single(interp1(tu,u,tu2,'spline'));
usig2 = single(interp1(tu,usig,tu2,'spline'));

check1 = abs(u2)>2500;
gate = check1;

check2 = shiftdata(check1,95,'circular'); % 75 is centered.

phase_shift = check1.*check2;

% edge knock off filter
% [num_f,den_f]=ChebyshevFilter(.05,0,0,4);
% [yout] = applyIIR(num_f,den_f,phase_shift);
% plot([phase_shift' yout']);
% return;


 coho    = single(exp(1i*2*pi*fif_1*tu2));
 coho_15    = single(exp(1i*2*pi*15.255e6*tu2));
 coho_4p397    = single(exp(1i*2*pi*4.397e6*tu2));
 sum_upc = coho.*usig2;
 rng_upc = coho.*usig2;
 rng_upc = rng_upc .* exp(1i*pi*(phase_shift));
 figure;
 subplot(2,1,1);
 plot(tu2,real(sum_upc));ah = gca;grid on;
 subplot(2,1,2);
 plot(tu2,real(rng_upc));ah = [ah gca];grid on;
 linkaxes(ah,'xy');
 %xlim([0 14e-7]);
 set(gcf,'Name','Sum&Rng');
 
 % 
 %
%  [bpf_124_num,bpf_124_den]=ChebyshevFilter( ...
%     (12.0954e6+.1237e6/2)/(1/dt2), ... % Break
%     0, ...      % 0 = Low Pass, 1 = High Pass
%     0, ...      % Percent Ripple (0 to 29)
%     16);         % Number of poles
% 
%  [cout,frq]=dbode(bpf_124_num,bpf_124_den,1,dt2,1000);
%  figure;
%  plot(1e-6*frq,10*log10(abs(cout)));
%  xlim([11.0 13]);
 
 flt=firFilter(12.03355e6,12.15725e6,1/dt2,256*32,'window_cheby(x,80)');
 sum_upc_flt = applyFilter(flt,sum_upc);
 rng_upc_flt = applyFilter(flt,rng_upc);
 
 sum_2 = sum_upc_flt ./ coho_15;
 rng_2 = rng_upc_flt ./ coho_15;
 
 %figure;
 %plot(real([sum_2.' rng_2.']))
 
  flt2=firFilter(3.2e6-0.5*1.61e3,3.2e6+0.5*1.61e3,1/dt2,256*256,'window_cheby(x,140)');
 sum_2_flt = applyFilter(flt2,sum_2);
 rng_2_flt = applyFilter(flt2,rng_2);
 
 figure;
 plot(real([sum_2_flt.' rng_2_flt.']));
 set(gcf,'Name','Sum&Rng 2nd IF');
 
 sum_3 = sum_2_flt .* coho_4p397; % multiply here because we are sum/rng are in "negative freq" space
 rng_3 = rng_2_flt .* coho_4p397;
 
 figure;
 plot(real([sum_3.' rng_3.']));
 set(gcf,'Name','Sum&Rng 3rd IF'); 
 
 sum_3_ddc = sum_3(1:10:end);
 rng_3_ddc = rng_3(1:10:end);
 dt3 = dt2*10;
 tu3 = tu2(1:10:end);
 
 flt3=firFilter(1.2e6-0.5*240,1.2e6+0.5*240,1/dt3,256*256,'window_cheby(x,140)');
 sum_3_ddc_flt = applyFilter(flt3,sum_3_ddc);
 rng_3_ddc_flt = applyFilter(flt3,rng_3_ddc);
 
 figure;
 plot(real([sum_3_ddc_flt.' rng_3_ddc_flt.']));
set(gcf,'Name','Sum&Rng 3rd IF BPF');

%%
filterBandwidth = 194.820e3;
filterCenter = 4*filterBandwidth;
filterHigh   = filterCenter + filterBandwidth/2;
filterLow    = filterCenter - filterBandwidth/2;

fracSampleRateBandwidth = filterHigh/((1/dt));
highpassFlag = 0;
[num_rough,den_rough]=ChebyshevFilter(fracSampleRateBandwidth,highpassFlag,0,6);
 [output,Fd] = dbode(num_rough,den_rough,1,dt,1000);
 %output = abs(output);
 figure;
 plot(Fd/1e6,20*log10(output));
 xlim(1e-6*[filterCenter-5*filterBandwidth filterCenter+5*filterBandwidth]);

 %

%function [num,den] = ChebyshevFilter(FC,LH,PR,NP)

 [yout] = applyIIR(num_rough,den_rough,u);
 figure;
 plot([abs(yout)' abs(u)]);

 xlim([81688.0990336257          81790.8298034776]);
 ylim([-3928.30875757902          23191.5742833567]);