close all;
clear all;
clc;
matlab_root = 'I:/MATLAB/';
addpath([matlab_root 'Analysis']);
addpath([matlab_root 'Antenna/Array_BLUE']);
addpath([matlab_root 'Printing']);
addpath([matlab_root 'Windows']);
addpath([matlab_root 'DSP']);
addpath([matlab_root 'DigitalFilters/Iowa Hills MATLAB']);
addpath([matlab_root 'Signal Generation/PWM Signal']);
addpath([matlab_root 'Signal Generation/Noise Signal']);
pp = PrepForPrint();
ON = 1;
OFF = 0;
TRUE = 1;
FALSE = 0;

if(0)
    %%
    for k = 1 : 6
        p = get(k,'Position');
        fprintf('set(%d,''Position'',[%f %f %f %f]);\n',k,p)
    end
end

% Instantiate array object
aesa = array('4-Element Array');
RF = 3e8/(2*sqrt(pi));
%RF = 10e9;
c = 3e8;
lambda = c / RF;
nzFigDb = 6;
nzFig = 10.^(nzFigDb / 10);
kboltz = 1.38e-23;

angleSteerDeg = [0 0]; % az x el

aesa.CreateRectangularArray( ...
    c/RF, ...   % wavelength (m)
    1.0,  ...   % efficiency (%/100)
    0.75, ...   % dx spacing x lambda
    0.75, ...   % dy spacing x lambda
    4,    ...   % # of x direction elements
    1);         % # of y direction elements

% set amplitude & phase noise
aesa.aVar = 0*0.5;
aesa.pVar = 0*10*pi/180;

elGainDbi = 0;
elGainDbi_V = elGainDbi / 2.0;
elGain_V = 10^(elGainDbi_V/10);
aesa.upd_elementGaindBV =elGainDbi_V;

aesa.recenterGeometry();
aesa.upd_redefineCoordinates();
aesa.upd_scaleWeights();

%% Generate simple dwell
Fs   = 1e6;
nzBW = Fs; 
Fif  = 0e6;
Pri  = 4e-3;
CPI  = 400e-3;
t    = single(0:(1/Fs):CPI);
N = length(t);

% create four correlated signal sources
nzAmp = 1.0e-0 / sqrt(2) * sqrt(kboltz*290*Fs*nzFig);
extNzAmp = 1.0e-4 / sqrt(2);
fprintf('Ext Nz Amp (dBW): %f\n',20*log10(extNzAmp)+3);
fprintf('Sub Nz Amp (dBW): %f\n',20*log10(extNzAmp)+3+12);
fprintf('Indiv. Channel Noise Floor (dBW): %f\n',20*log10(nzAmp)+3);
fprintf('Sum    Channel Noise Floor (dBW): %f\n',20*log10(nzAmp*sqrt(4))+3); % 4 elements, +3 for /sqrt(2) above
fprintf('Sum Channel SNR (dB): %f\n',20*log10(extNzAmp/nzAmp)+6);
%s1 = extNzAmp * (randn(1,N) + 1i*randn(1,N));
s1 = (0.1*exp(1i*2*pi*50e3*t));
s2 = s1;
s3 = s1;
s4 = s1;

% create four uncorrelated noise sources
nz1 = nzAmp * (randn(1,N) + 1i*randn(1,N));
nz2 = nzAmp * (randn(1,N) + 1i*randn(1,N));
nz3 = nzAmp * (randn(1,N) + 1i*randn(1,N));
nz4 = nzAmp * (randn(1,N) + 1i*randn(1,N));

% TOTAL signal voltage gain through summing junction is:
% (1.0 V/m) x (net unity) x (# elements) = 4.0 V (assuming R = 1 Ohm)
% TOTAL noise voltage gain is Vn,rms x sqrt(# elements) = 

% Verify
sumJunc_nz = nz1+nz2+nz3+nz4;
sumJunc_sig = s1+s2+s3+s4;
sumJunc = sumJunc_sig+sumJunc_nz;

wo_nz1 = WaveObj(nz1,[],Fs);
wo_s1 = WaveObj(s1,[],Fs);

wo_nz = WaveObj(sumJunc_nz,[],Fs);
wo_sig = WaveObj(sumJunc_sig,[],Fs);
wo_sum = WaveObj(sumJunc,[],Fs);
%%
plot(wo_nz1);
subplot(2,1,1);
title('Single Channel Noise PSD');
add_print_callbacks;

plot(wo_s1);
subplot(2,1,1);
title('Single Channel Signal PSD');
add_print_callbacks;

plot(wo_sig);
subplot(2,1,1);
title('Summation Signal PSD');
add_print_callbacks;

plot(wo_nz);
subplot(2,1,1);
title('Summation Noise PSD');
add_print_callbacks;

% Signal is increased by 16x (power) | 4x (voltage)
% Noise increased by 4x (power) | 2x (voltage)
% Signal / Noise in Voltage is improved by 2x
% getGain gives a value of 2.0

%% Now check subapertues

% clear sub-array weights
aesa.upd_setSubArrayWeight(1);

%initialize/zero
subArray_iq = 0;


nSubArrays = length(unique(aesa.upd_subArrayNumber));

% add this source signal to elemental sampling
subArray_iq = subArray_iq + ...
    aesa.upd_genIQ_Channels( ...
    s1, ...
    'fromAz', 0, ...
    'fromEl', 0, ...
    'subArray', 1:nSubArrays, ...
    'applySteering', 1, ...
    'sumOutput', 0);

subArray_iq_nz = 0*subArray_iq;

% Source signal generation / propagation complete - now add rx noise + FIR
% filter to match bandwith
for ks = 1 : nSubArrays
    subArray_iq_nz(ks,:) = subArray_iq_nz(ks,:) + nzAmp * (randn(1,N) + 1i*randn(1,N));
end

sumJunc2_nz = sum(subArray_iq_nz,2);
sumJunc2_sig = sum(subArray_iq,2);
sumJunc = sumJunc2_sig+sumJunc2_nz;

wo_nz2 = WaveObj(sumJunc_nz,[],Fs);
wo_sig2 = WaveObj(sumJunc_sig,[],Fs);
wo_sum2 = WaveObj(sumJunc,[],Fs);

plot(wo_sig2);
subplot(2,1,1);
title('SA Summation Signal PSD');
add_print_callbacks;

plot(wo_nz2);
subplot(2,1,1);
title('SA Summation Noise PSD');
add_print_callbacks;

% Same result.

%%


%%
Rest = aesa.upd_estCorrMtx(1*subArray_iq+1*subArray_iq_nz,[]);
            
            %
[Vc,Dc]=eig(Rest,'vector');
[~,sortIdx]=sort(Dc,'descend');
Dc = Dc(sortIdx);
Vc = Vc(:,sortIdx);

10*log10(abs(Dc))

set(1,'Position',[236.000000 66.000000 560.000000 420.000000]);
set(2,'Position',[235.000000 596.000000 560.000000 420.000000]);
set(3,'Position',[799.000000 596.000000 560.000000 420.000000]);
set(4,'Position',[798.000000 66.000000 560.000000 420.000000]);
set(5,'Position',[1361.000000 594.000000 560.000000 420.000000]);
set(6,'Position',[1360.000000 65.000000 560.000000 420.000000]);
