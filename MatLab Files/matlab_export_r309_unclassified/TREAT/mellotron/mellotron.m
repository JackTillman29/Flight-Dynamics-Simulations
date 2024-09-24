close all;
clear;
clc;
addpath('I:\MATLAB\TREAT\test_tools\scope_data_view_score_gui');
matlab_root = 'I:\MATLAB';
addpath([matlab_root '\Analysis']);
addpath([matlab_root '\Printing']);
addpath([matlab_root '\Visualization']);

pp = PrepForPrint();

dataFolder = '20181205_152459_Lancelot\Trial0_0';
%dataFolder = '20181205_150357_Mel_235_15_BPOFF\Trial0_0';
%dataFolder = '20181205_151553_Mel_235_15_BPON\Trial0_0';
%dataFolder = '20181205_151750_Mel_50_200_BPON\Trial0_0';
%dataFolder = '20181205_150110_Mel_50_200_BPOFF\Trial0_0';

inChanStr  = 'scope_channel_2.bin';
outChanStr = 'scope_channel_4.bin';

inData = oscope_channel_string_and_binary_reader([dataFolder '\' inChanStr]);
outData = oscope_channel_string_and_binary_reader([dataFolder '\'  outChanStr]);
inData.t = inData.t - inData.t(1);

wo_in = WaveObj(inData.Vc,[],inData.Fs,length(inData.t));
wo_out = WaveObj(outData.Vc,[],inData.Fs,length(inData.t));


%% Pulse data processing
% Find when input pulse envelope exceeds threshold
pulsePresThresh = max(abs(inData.Vc)) / 2.0; % can set manually if needed
inSigAboveThresh = abs(inData.Vc) > pulsePresThresh;

outputPresThresh = max(abs(outData.Vc)) / 2.0; % can set manually if needed
outSigAboveThresh = abs(outData.Vc) > outputPresThresh;

idxPulseStart = find(diff(inSigAboveThresh) == 1);
nPulseStart = length(idxPulseStart);
fprintf('Estimated # of Pulses: %d\n',nPulseStart);

priPerLine = mean(diff(idxPulseStart)) * inData.dt;
prfEst = 1./priPerLine; % for later
fprintf('Estimated PRI (us): %18.6f\n',1e6*priPerLine);
priPerLine = 1.05 * priPerLine;
smpPerLine = floor(priPerLine ./ inData.dt);

%% Coho ID and tuning
[maxf,maxfi]=max(abs(wo_in.fft)); % get a "close" guess from fft
cohoEst = wo_in.frq(maxfi); % based on max of input fft
fprintf('Estimated Coho (FFT): %18.6f\n',cohoEst);
cohoData = single(exp(1i*(2*pi*cohoEst)*inData.t));

cohoTestData = angle(inData.Vc./cohoData); 
cohoTestData(inSigAboveThresh == 0) = nan; % only data where input is present

% Pass "pulse present" phase data to polyfit to compute residual
fineTune = polyfit( ...
    inData.t(~isnan(cohoTestData)), ...
    cohoTestData(~isnan(cohoTestData)), ...
    1); % First order (assumes no LO drift)

fprintf('Estimated Coho (PF): %18.6f\n',cohoEst+2*pi*fineTune(1));
cohoData = single(exp(1i*(2*pi*cohoEst+fineTune(1))*inData.t));

%% 1-D Plots
dataFolderClean = strrep(dataFolder,'\','\\');
dataFolderClean = strrep(dataFolderClean,'_','\_');

figure;

plot(inData.t(outSigAboveThresh ~= 0), ...
    angle(outData.Vc(outSigAboveThresh ~= 0)./cohoData(outSigAboveThresh ~= 0)),'.');
hold on;
plot(inData.t(inSigAboveThresh ~= 0), ...
    angle(inData.Vc(inSigAboveThresh ~= 0)./cohoData(inSigAboveThresh ~= 0)),'.');
hold off;
xlabel('Time (s)');
ylabel('rad');
title({'Coherence Check',['\fontsize{8}' dataFolderClean]});
grid on;
legend('Output','Input');
add_analysis_callbacks;
PrepForPrint(gcf,pp);add_print_callbacks;
set(gcf,'Position',[1131.7          437          560          420]);


figure;
hs=subplot(2,1,1);
plot(hs(1),1e6*inData.t,real(inData.Vc),'Color',0.8*[1 1 1]);
set(hs(1),'NextPlot','add');
plot(hs(1),1e6*inData.t,[abs(inData.Vc)]);
xlabel('Time (\mus)');ylabel('mV');
grid on;title('Input');
hs(2) = subplot(2,1,2);
plot(hs(2),1e6*inData.t,real(outData.Vc),'Color',0.8*[1 1 1]);
set(hs(2),'NextPlot','add');
plot(hs(2),1e6*inData.t,[abs(outData.Vc)]);
linkaxes(hs,'x');
xlabel('Time (\mus)');ylabel('mV');
grid on;
title('Output');
PrepForPrint(gcf,pp);add_print_callbacks;
set(gcf,'Position',[1138.3       40.111          560       354.22]);

plot([wo_in wo_out],1e-3,'kHz',-1e-3*(cohoEst+fineTune(1)/(2*pi)))
PrepForPrint(gcf,pp);
add_print_callbacks;

%% 2-D Plots


% repackage data
Min = complex(zeros(nPulseStart,smpPerLine,'single'));
Mout = complex(zeros(nPulseStart,smpPerLine,'single'));
Mcoho = complex(zeros(nPulseStart,smpPerLine,'single'));
Mpp = zeros(nPulseStart,smpPerLine,'single');

for k = 1 : nPulseStart
    try
        Min(k,:) = inData.Vc(idxPulseStart(k) : idxPulseStart(k) - 1 + smpPerLine);
        Mout(k,:) = outData.Vc(idxPulseStart(k) : idxPulseStart(k) - 1 + smpPerLine);
        Mcoho(k,:) = cohoData(idxPulseStart(k) : idxPulseStart(k) - 1 + smpPerLine);
    catch
        if(k == nPulseStart)
            n = length(inData.Vc(idxPulseStart(k) : end));
            Min(k,1:n) = inData.Vc(idxPulseStart(k) : end);
            Mout(k,1:n) = outData.Vc(idxPulseStart(k) : end);
            Mcoho(k,1:n) = cohoData(idxPulseStart(k) : end);
        else
            warning(['Data remap error on pulse ' num2str(k) ' of ' num2str(nPulseStart)]);
        end
    end
end



Mpp = abs(Min) > pulsePresThresh; % 1 or 0 if above threshold. Useful as a mask
%%
fastTimeVec = (0:size(Min,2)-1) * inData.dt * 1e6;
slowTimeVec = (inData.t(idxPulseStart) - inData.t(1)) * 1e6;

figure;
imagesc(fastTimeVec,slowTimeVec,angle(Min./Mcoho));
xlabel('Fast Time (\mus)');ylabel('Slow Time (\mus)');
title({'Input Coherent Phase Map',['\fontsize{8}' dataFolderClean]});
colorbar;ah =gca;
caxis([-pi pi]);PrepForPrint(gcf,pp);add_print_callbacks;
set(gcf,'Position',[5.4444  461.4444  560.0000  420.0000]);

figure;
imagesc(fastTimeVec,slowTimeVec,angle(Mout./Mcoho));
xlabel('Fast Time (\mus)');ylabel('Slow Time (\mus)');
title({'Output Coherent Phase Map',['\fontsize{8}' dataFolderClean]});
colorbar;ah = [ah gca];
caxis([-pi pi]);PrepForPrint(gcf,pp);add_print_callbacks;
set(gcf,'Position',[568.1111  461.4444  560.0000  420.0000]);

figure;
imagesc(fastTimeVec,slowTimeVec,1e3*abs(Min./Mcoho));
xlabel('Fast Time (\mus)');ylabel('Slow Time (\mus)');
title({'Input Amplitude Map (mV)',['\fontsize{8}' dataFolderClean]});
colorbar;ca = caxis();ah = [ah gca];
PrepForPrint(gcf,pp);add_print_callbacks;
set(gcf,'Position',[14.3333   49.0000  560.0000  339.5556]);

figure;
imagesc(fastTimeVec,slowTimeVec,1e3*abs(Mout./Mcoho));
xlabel('Fast Time (\mus)');ylabel('Slow Time (\mus)');
title({'Output Amplitude Map (mV)',['\fontsize{8}' dataFolderClean]});
colorbar;caxis(ca);ah = [ah gca];
PrepForPrint(gcf,pp);add_print_callbacks;
set(gcf,'Position',[576.1111   48.5556  560.0000  336.4444]);



%% RDMs
slowFreqVec = slowTimeVec ./ max(slowTimeVec) * prfEst;
slowFreqVec = slowFreqVec - prfEst / 2;


q = Min./Mcoho;
q(isnan(q)) = 0;
figure;
imagesc(1e-3*slowFreqVec,fastTimeVec,20*log10(abs(fftshift(fft(q,[],1)./size(q,1)))).');
xlabel('Frequency (kHz)');
ylabel('Time (\mus)');
title('Fast Time Doppler Map [Input]');
colorbar;colormap(colormap_fade2black);ca = caxis();
PrepForPrint(gcf,pp);add_print_callbacks;ah = gca;

q = Mout./Mcoho;
q(isnan(q)) = 0;
figure;
imagesc(1e-3*slowFreqVec,fastTimeVec,20*log10(abs(fftshift(fft(q,[],1)./size(q,1)))).');
xlabel('Frequency (kHz)');
ylabel('Time (\mus)');
title('Fast Time Doppler Map [Output]');
colorbar;colormap(colormap_fade2black);caxis(ca);
PrepForPrint(gcf,pp);add_print_callbacks;ah = [ah gca];
linkaxes(ah,'xy');
