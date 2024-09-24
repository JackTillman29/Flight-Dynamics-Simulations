close all;
clear all;
clc;

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
nSamples = 1000;

noise_azimuths = [-30 20]*pi/180;
noise_elevations = [0 0]*pi/180;
element_iq = 0;

aesa.upd_setSteer(0*pi/180,0);

for k = 1 : length(noise_azimuths)
    noise_sig = randn(1,nSamples) + 1i * randn(1,nSamples);
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


JNR_dB = 10; % dB
JNR = 10^(-JNR_dB/20);
aesa.pVar = 0*5*pi/180;

% add some temporal noise
element_iq = element_iq + JNR .* (randn(size(element_iq)) + 1i*randn(size(element_iq)));

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

if(0) % old method
    

    Dinvc = aesa.diagInv(Dc,length(noise_azimuths));
    Dinvc(Dinvc ~= 0) = 1; % each usable eigenvector will project as a unit vector

    % Optimal weighting vector is R^-1 * s
    %w = (V * Uinv * V') * s;
    
    %s = ones(aesa.nElements,1);
    wc = (Vc * Dinvc * Vc') * (s);

    %abf_termc = (s-wc)';
    abf_termc = conj(s-wc);

else % keith's new understanding
    sh = conj(s); % retro-directivity vector (phase front of incident wave)
    
    % Interference subspace
    V = Vc(:,1:length(noise_azimuths));
    
    % Find nearest steering vector conjugate that lies in the nullspace of
    % V
    shprime = (eye(aesa.nElements) - V*V') * sh;
    
    abf_termc = conj(shprime);
    
end
aesa.upd_plot( ...
    linspace(-pi/2,pi/2,2000), ...
    0, ...
    aesa.ENUM_MODE_UNIFORM);
YL = ylim();
ylim(YL);
aesa.abf_term = abf_termc;
aesa.upd_setSubArrayWeight(abf_termc);        

aesa.upd_plot( ...
    linspace(-pi/2,pi/2,2000), ...
    0, ...
    aesa.ENUM_MODE_UNIFORM);

hold on;
ylim(YL);
line( ...
    [noise_azimuths;noise_azimuths]*180/pi, ...
    repmat(YL',1,length(noise_azimuths)), 'Color','r','LineStyle','--');


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

aesa.upd_setSteer(testAngleAz,testAngleEl);

% Steering Vector (set via upd_setSteer)
sh = conj(exp(1i*aesa.phs_steer));

% What is this doing? 
% 1 - Project the steering vector (sh) onto the signal subspace
% 2 - Take that resultant vector in the signal subspace and project it onto
%     the steering vector.
% 
% The signal subspace spans the same space as steering vectors to all
% external interference source. (e.g. steering vectors toward jammers are
% linearly dependent on signal subspace eigenvectors). If a steering vector
% corresponds to an external jammer, it will project fully onto the signal
% subspace with zero loss in projection. 
%
% What we are interested in is only that component of the steering vector
% which lies in the signal subspace, then projected back onto the steering 
% vector, s. This informs how much of our steering vector lies within the 
% signal subspace.
signalSubspaceProjection = sh'*(Vs*Vs')*sh;

% By nature of the eigenvector decomposition of a positive definite hermitian
% matrix, all eigenvectors are orthonormal. This means they are all unit
% length and all orthogonal to eachother. What follows is that a similar 
% analysis can be done above using the noise subspace. In this case,
% steering vectors which project into the null space of the noise subspace
% must be entirely representable in the signal subspace. In other words,
% steering vectors which project to vector zero in the noise subspace will
% correspond to signal direction of arrival.

% What is the noise subspace? Its the space spanned by the orthonormal 
% eigenvector set associated with the small (noise-like) eigenvalues.  
% Remember that a full span (Signal & Noise) can be represented as a vector
% sum of two vectors: 1 - component in the noise subspace + 2 - component in
% signal subspace.
%
% e.g. signal = (proj of signal on noise subspace [Vn]) + ...
%               (proj of signal on signal subspace [Vs])
%
%  S = Pn(S) + Ps(S) = Vn*Vn'*S + Vs*Vs'*S
%  
% The noise subspace projection is
% Vn*Vn'S = S*(I - Vs*Vs')
%
% This can be found in Farina pp. 315 for reference.
noiseSubspaceProjection1 = sh'*(Vn*Vn')*sh; % Method 1
noiseSubspaceProjection2 = sh'*(eye(aesa.nElements)-Vs*Vs')*sh; % Method 2

% should be large because we are pointed right at a source
MUSIC = 1./abs(noiseSubspaceProjection1)

%sa_music_wgt.'*(Vc(:,(nSigEig+1):end)*Vc(:,(nSigEig+1):end)')*sa_music_wgt

%% Compute MUSIC for array of far-field
nAz = 400;
nEl = 400;
azVec = linspace(-60,60,nAz);
elVec = linspace(-60,60,nEl);
% sh will be [N elements x (N Az x N El)]
sha = conj(exp(1i*aesa.upd_setSteerArray( ...
    azVec*pi/180, ...
    elVec*pi/180)));

% project all steering vectors onto noise subspace (method 2)
% NOTE: THis is done slightly differently than above because x'*x does not
% yeild a scalar dot product for the full matrix implementation. This
% requires the dot product to be implemnted as an element-wise multiply and
% summation. They are equivalent.
tmp = (eye(aesa.nElements)-Vs*Vs')*sha;
noiseSubProjArr2 = reshape(sum(conj(sha).*tmp,1),nEl,nAz);

MUSIC_Arr = 10*log10(1./abs(noiseSubProjArr2));

figure;
imagesc(azVec,elVec,MUSIC_Arr);
xlabel('Azimuth Angle of Arrival');
ylabel('Elevation Angle of Arrival');
title('MUSIC Az/El Response [dB]');
colorbar;
set(gca,'YDir','normal');

% 1D
nAz = 10000;
nEl = 1;
azVec = linspace(-90,90,nAz);
sh1D = conj(exp(1i*aesa.upd_setSteerArray( ...
    azVec * pi/180, ...
    0)));

% project all steering vectors onto noise subspace (method 2)
tmp = (eye(aesa.nElements)-Vs*Vs')*sh1D;
noiseSubProj1D2 = reshape(sum(conj(sh1D).*tmp,1),nEl,nAz);

MUSIC_1D = 1*(1./abs(noiseSubProj1D2));
figure;
plot(azVec,MUSIC_1D,'.-');
grid on;
xlabel('Angle of Arrival');
ylabel('Linear');
title('MUSIC Response');
