close all;
clear all;
clc;
matlab_root = 'I:/MATLAB/';
addpath([matlab_root 'Analysis']);
addpath([matlab_root 'Printing']);
addpath([matlab_root 'Windows']);
addpath([matlab_root 'DigitalFilters/Iowa Hills MATLAB']);
addpath([matlab_root 'Signal Generation/PWM Signal']);
addpath([matlab_root 'Signal Generation/Noise Signal']);
pp = PrepForPrint();

% Notes:
%  N = # of elements

% Instantiate array object
aesa = array('Sample Array');

% Initialize the array object
aesa.InitStapArray( ...
    3e8/34e9, ...               % wavelength, in meters
    1.0, ...                    % efficiency, between 0.0 and 1.0
    'Sample_Elements.txt', ...  % element position file, x,y,z meters [3 x N]
    [], ...                     % amplitude sum taper file, [N x 1]
    [], ...                     % amplitude delta az taper file, [N x 1]
    [], ...                     % amplitude delta el taper file, [N x 1]
    []);                        % element pattern file, not implemented in InitStapArray

% Create the array taper using internal routines
aesa.compTaylorWgt(4,-30);      % compute taylor (rx) weights
aesa.compBaylissWgt(4,-40);     % compute difference pattern weights

% set amplitude & phase noise
%aesa.aVar = 0*0.5;
aesa.pVar = 1*10*pi/180;

% set element pattern (NOTE: This simply creates a 1-D solid angle cos^2
% element pattern with a peak power gain of elGainDbi).
elGainDbi = 4.5;
elGainDbi_V = elGainDbi / 2.0;
elGain_V = 10^(elGainDbi_V/10);
aesa.matlab_setElementPattern( ...
    linspace(-pi/2,pi/2,100), ...
    elGain_V .* cos(linspace(-pi/2,pi/2,100)));


aesa.upd_redefineCoordinates();
aesa.upd_scaleWeights();


%% Create sub-array definition
%aesa.upd_subArrayNumber = ...
%    [ones(1,32 )*1 ones(1,96)*2 ones(1,96)*3 ones(1,128)*4 ...
%     ones(1,128)*5 ones(1,96)*6 ones(1,96)*7 ones(1,32 )*8 ]';

% sub array
load new_sub_original.mat;
for ksub = 1 : 16
    %    yd = eval(sprintf('new_sub%d(:,2);',ksub));
    subArrayIdx = eval(sprintf('new_sub%d(:,1);',ksub));
    % hp=plot(yd,zd,'.','MarkerSize',18);
    % disp(length(zd))
    aesa.upd_subArrayNumber(subArrayIdx) = ksub;
end
clear subArrayIdx new_sub*;

aesa.upd_setSubArrayWeight(ones(16,1));
aesa.plotGeometry(aesa.upd_subArrayNumber);

%% Generate simple PD dwell
Fs = 25e6;
Fif = 0e6;
Pd = 0.75e-6;
Prf = 600e3;
Pri = 1 / Prf;
CPI = 4e-3;
t = single(0:(1/Fs):CPI);
Pin = 0; % dBm
Ain = 10.^((Pin-30)./20);

%Psnr_element_dB = -20%-20-33;
%Psnr_element = 10.^(Psnr_element_dB./10);

% Generate simple unmodulated pulsed Doppler waveform
sig_in = Ain*exp(1i*2*pi*Fif*t) .* PWM_Signal( ...
    Pd, ...
    Pri, ...
    1/Fs, ...
    length(t), ...
    1);

% PRF line filter
ih = FirFilterC();
ih.AutogenFir(Fs,2048,-10e3,10e3);

% Truncate data such that it all fits nicely in 2D FTST map
nPri = round(Pri * Fs);
nCPI = floor(length(sig_in) / nPri);
sig_in = sig_in(1:(nCPI*nPri));
t = t(1:(nCPI*nPri));



%% create a noise source out there
% now we use the PD data as if it is coming in on multiple array channels


nSamples = 100000;
nSubArrays = 16;


jam_az = 0;
angleSteerDeg = [0 0];
for kjam = 1 : length(jam_az)
    noise_azimuths = [jam_az(kjam)]*pi/180;
    noise_elevations = [0]*pi/180;
    
    
    
    % Apply full array phase shifting
    aesa.upd_setSteer(angleSteerDeg(1)*pi/180,angleSteerDeg(2)*pi/180);
    %if(kjam == 1)
    subArray_iq = 0;
    for k = 1 : length(noise_azimuths)
        %noise_sig = randn(1,nSamples) + 1i * randn(1,nSamples);
        %noise_sig = 1 + 1i*0;
        noise_sig = sig_in;
        noise_azimuth = noise_azimuths(k);
        noise_elevation = noise_elevations(k);
        
        % clear sub-array weights
        aesa.upd_setSubArrayWeight(1);
        
        subArray_iq = subArray_iq + ...
            aesa.upd_genIQ_Channels( ...
            noise_sig, ...
            'fromAz', noise_azimuth, ...
            'fromEl', noise_elevation, ...
            'subArray', 1:nSubArrays, ...
            'applySteering', 1, ...
            'sumOutput', 0);
    end
    %end
    
    % now add noise to sub-array iq data
    for ks = 1 : nSubArrays
        subArray_iq(ks,:) = subArray_iq(ks,:) + ...
            1*GenerateUnfilteredNoise(nPri*nCPI,Fs,1e-4,'gaussian_iq');
        
    end
    
    % now sum all data to form Sum channel (for reference)
    sumChannel = sum(subArray_iq,1);
    
    %%
    FTST = reshape(sumChannel,nPri,nCPI).';

    % Compute RDM
    RDM = fftshift(fft(FTST,[],1)./nCPI).';
    fvec = 1e-3*linspace(-Prf/2,Prf/2,nCPI);
    rvec = 150e6*linspace(-Pri/2,Pri/2,nPri);

    fig_rdm = figure;
    imagesc(fvec,rvec,20*log10(abs(RDM)));
    xlabel('Doppler [kHz]');
    ylabel('Range [m]');
    set(gca,'YDir','normal');
    title({'Array Range Doppler Matrix',sprintf('PRF=%3.1fkHz | PD=%3.2f\\mus',Prf/1000,Pd*1e6)});
    xlim([-10 10]);
    add_analysis_callbacks;
    
    %%
    % simulate filtering SA data to rid PRF lines
    if(1)
    for ks = 1 : nSubArrays
        subArray_iq(ks,:) = ih.applyFirData(subArray_iq(ks,:));
    end
    end
    
    
    % Estimate correlation matrix
    Rest = aesa.upd_estCorrMtx(subArray_iq,[]);
    
    %
    [Vc,Dc]=eig(Rest);
    [~,sortIdx]=sort(Dc,'descend');
    Dc = diag(Dc(sortIdx));
    Vc = Vc(:,sortIdx);
    %eigIdx = length(noise_azimuths);
    eigIdx =sum(abs(Dc) > 1e-6);
    V_is = Vc(:,1:eigIdx);
 
   
    % Now perform MUSIC processing
    % Compute MUSIC for array of far-field rel main beam
    nAz = 400;
    nEl = 400;
    azVec = linspace(-5,5,nAz);
    elVec = linspace(-5,5,nEl);
    
    % sh (sha) will be [N subarrays x (N Az x N El)]
    sha = conj( exp(1i*aesa.upd_setSteerArray( ...
                azVec*pi/180, ...
                elVec*pi/180, ...
                'subarray')));
            
    % project all steering vectors onto noise subspace (method 2)
    % NOTE: THis is done slightly differently than above because x'*x does not
    % yeild a scalar dot product for the full matrix implementation. This
    % requires the dot product to be implemnted as an element-wise multiply and
    % summation. They are equivalent.
    tmp = (eye(nSubArrays)-V_is*V_is')*sha;
    
    noiseSubProjArr2 = reshape(sum(conj(sha).*tmp,1),nEl,nAz);
    
    MUSIC_Arr = 10*log10(1./abs(noiseSubProjArr2));
    if(kjam == 1)
        figure;
        him=imagesc(azVec,elVec,MUSIC_Arr);
        xlabel('\DeltaAzimuth Angle of Arrival');
        ylabel('\DeltaElevation Angle of Arrival');
        title('MUSIC Az/El Response [dB]');
        colorbar;
        set(gca,'YDir','normal');
        set(gcf,'Position',[1137          48         560         420]);
        %caxis([-10 10]);
    else
        set(him,'CData',MUSIC_Arr);
    end
    
%     return;
end