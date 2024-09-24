close all;
clear all;
clc;

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
aesa.aVar = 0*0.5;
aesa.pVar = 0*10*pi/180;

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



%aesa.upd_setPhaseQuantBits(5);

aesa.upd_setSteer(-20*pi/180,0*pi/180);

aesa.upd_plot( ...
    linspace(-pi/2,pi/2,1000), ...
    0, ...
    aesa.ENUM_MODE_UNIFORM);


%% create a noise source out there
nSamples = 1000;

noise_azimuths = [18 30]*pi/180;
noise_elevations = [0 0]*pi/180;
element_iq = 0;

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



% Estimate correlation matrix
Rest = aesa.upd_estCorrMtx(element_iq,[]);


%
[Vc,Dc]=eig(Rest);
[~,sortIdx]=sort(Dc,'descend');
Dc = diag(Dc(sortIdx));
Vc = Vc(:,sortIdx);
eigIdx = length(noise_azimuths);
V_is = Vc(:,1:eigIdx);
%%
sh = exp(1i*aesa.phs_steer*0); % steering applied already (this is relative)
shprime = (eye(aesa.nElements)-V_is*V_is') * sh;
abf_termc = conj(shprime);
aesa.upd_setSubArrayWeight(abf_termc);  


% Plot antennas
%aesa.upd_plot(linspace(-pi/2,pi/2,5000),0, ...
    %[aesa.ENUM_MODE_UNIFORM aesa.ENUM_MODE_ABF],{'Uniform','Adaptive'});
    
aesa.upd_plot( ...
    linspace(-pi/2,pi/2,300), ...
    linspace(-pi/2,pi/2,300), ..., ...
    aesa.ENUM_MODE_UNIFORM);