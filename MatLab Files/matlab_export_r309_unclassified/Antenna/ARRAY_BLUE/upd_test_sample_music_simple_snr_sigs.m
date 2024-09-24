close all;
clear all;
clc;
matlab_root = 'I:\MATLAB\';
addpath([matlab_root 'Signal Generation\LFM Signal']);
% Notes: 
%  N = # of elements

% Instantiate array object
aesa = array('5 Element ULA');

% Initialize the array object
aesa.CreateRectangularArray(3e8/34e9, 1.0, .5, .5, 5, 5);

% set element pattern (NOTE: This simply creates a 1-D solid angle cos^2
% element pattern with a peak power gain of elGainDbi).
elGainDbi = 0;
elGainDbi_V = elGainDbi / 2.0;
elGain_V = 10^(elGainDbi_V/10);
aesa.matlab_setElementPattern( ...
    linspace(-pi/2,pi/2,1000), ...
    elGain_V .* cos(linspace(-pi/2,pi/2,1000)));
aesa.upd_elementGaindBV = 0;

aesa.recenterGeometry();
aesa.upd_redefineCoordinates();

%aesa.upd_scaleWeights();

%aesa.upd_setSteer(-20*pi/180,0*pi/180);



%% create a noise source out there

nSamplesArray = logspace(1,6,6);
nSamplesArray = 100000;
%nSamplesArray = 1e7;
snrArray = linspace(-40,-10,100);
%snrArray = linspace(-60,-30,30);

res.peak = 0*nSamplesArray;
res.bottom = 0*nSamplesArray;
res.nSamples = 0*nSamplesArray;
res.snr = 0*nSamplesArray;
for ktrade_smp = 1 : length(nSamplesArray)
for ktrade_snr = 1 : length(snrArray)
nSamples = nSamplesArray(ktrade_smp)


noise_azimuths = [0]*pi/180;
noise_elevations = [0 0]*pi/180;
element_iq = 0;

aesa.upd_setSteer(0*pi/180,0);

JNR_dB = snrArray(ktrade_snr); % dB at the element level
JNR = 10^(-JNR_dB/20);
aesa.pVar = 0*5*pi/180;

[y,t,freq] = SimpleLfm(-2.5e6,2.5e6,128e-6,128e-6/nSamples);

for k = 1 : length(noise_azimuths)
    noise_sig = randn(1,nSamples) + 1i * randn(1,nSamples);
    %noise_sig = 1*exp(1i*2*pi./100 * (1:nSamples));
    %noise_sig = y;
    noise_azimuth = noise_azimuths(k);
    noise_elevation = noise_elevations(k);
    % Generate per-element I/Q data [Ne x #TD Samples]
    %element_iq = element_iq + aesa.upd_genIQ_Channels(noise_sig,noise_azimuth,noise_elevation);
    element_iq = element_iq + ...
           aesa.upd_genIQ_Channels( ...
            noise_sig, ...
            'fromAz', noise_azimuth, ...
            'fromEl', noise_elevation, ...
            'applySteering', 1, ...
            'sumOutput', 0);
end

% full gain SNR will be:
fullGainSNR = JNR_dB + 10*log10(aesa.nElements);
fprintf('Full Gain SNR (dB): %f\n',fullGainSNR);


% add some temporal noise
element_iq = element_iq + JNR .* sqrt(2)/2*(randn(size(element_iq)) + 1i*randn(size(element_iq)));

% Estimate correlation matrix
Rest = aesa.upd_estCorrMtx(element_iq,[]);


%
[Vc,Dc]=eig(Rest);
% Set are orthonormal 239 farina
[~,sortidxc]=sort(abs(diag(Dc)),'descend');
Dc = Dc(sortidxc,sortidxc);
Vc = Vc(:,sortidxc);
aesa.amp = ones(aesa.nElements,1);
s = exp(1i*aesa.phs_steer) .* aesa.amp;
sh = conj(s); % retro-directivity vector (phase front of incident wave)
disp([abs(Dc(1,1)) abs(Dc(2,2))])
% Interference subspace
V = Vc(:,1:length(noise_azimuths));
    

%% Let the music flow...
aesa.upd_setSubArrayWeight(1);

testAngleAz =  noise_azimuths(1);
testAngleEl =  noise_elevations(1);

% Real world: some algorithm should interpret the eigenvalues to identify
% the signal and noise subspace breakpoints.
nSigEig   = length(noise_azimuths);
nNoiseEig = length(Vc(:,1)) - nSigEig;

Vs = Vc(:,1:nSigEig);
Vn = Vc(:,nSigEig+1 : end);

%% Compute MUSIC for array of far-field
% nAz = 1;
% nEl = 400;
% azVec = linspace(-60,60,nAz);
% elVec = linspace(-60,60,nEl);
% % sh will be [N elements x (N Az x N El)]
% sha = conj(exp(1i*aesa.upd_setSteerArray( ...
%     azVec*pi/180, ...
%     elVec*pi/180)));

% project all steering vectors onto noise subspace (method 2)
% NOTE: THis is done slightly differently than above because x'*x does not
% yeild a scalar dot product for the full matrix implementation. This
% requires the dot product to be implemnted as an element-wise multiply and
% summation. They are equivalent.
% tmp = (eye(aesa.nElements)-Vs*Vs')*sha;
% noiseSubProjArr2 = reshape(sum(conj(sha).*tmp,1),nEl,nAz);
% 
% MUSIC_Arr = 10*log10(1./abs(noiseSubProjArr2));
% 
% figure;
% imagesc(azVec,elVec,MUSIC_Arr);
% xlabel('Azimuth Angle of Arrival');
% ylabel('Elevation Angle of Arrival');
% title('MUSIC Az/El Response [dB]');
% colorbar;
% set(gca,'YDir','normal');

% 1D
nAz = 1000;
nEl = 1;
azVec = linspace(-90,90,nAz);
sh1D = conj(exp(1i*aesa.upd_setSteerArray( ...
    azVec * pi/180, ...
    0)));

% project all steering vectors onto noise subspace (method 2)
tmp = (eye(aesa.nElements)-Vs*Vs')*sh1D;
noiseSubProj1D2 = reshape(sum(conj(sh1D).*tmp,1),nEl,nAz);

MUSIC_1D = 1*(1./abs(noiseSubProj1D2));
% figure;
% plot(azVec,10*log10(MUSIC_1D),'.-');
% grid on;
% xlabel('Angle of Arrival');
% ylabel('10Log');
% title('MUSIC Response');

res.peak(ktrade_snr) = max(MUSIC_1D);
res.bottom(ktrade_snr) = min(MUSIC_1D);
res.nSamples(ktrade_smp) = nSamples;
res.snr(ktrade_snr) = JNR_dB;



end
figure(1);
hold on;
plot(res.snr+13.9794,10*log10(res.peak./res.bottom));
drawnow;

figure(2);
hold on;
plot(res.snr,10*log10(res.peak./res.bottom));
drawnow;


end