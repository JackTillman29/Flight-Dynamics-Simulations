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
aesa.pVar = 0*5*pi/180;

% set element pattern (NOTE: This simply creates a 1-D solid angle cos^2
% element pattern with a peak power gain of elGainDbi).
elGainDbi = 4.5;
elGainDbi_V = elGainDbi / 2.0;
elGain_V = 10^(elGainDbi_V/10);
aesa.matlab_setElementPattern( ...
    linspace(-pi/2,pi/2,100), ...
    elGain_V .* cos(linspace(-pi/2,pi/2,100)));


aesa.upd_redefineCoordinates();

aesa.upd_setSteer(0*pi/180,10*pi/180);

aesa.upd_plot( ...
    linspace(-pi/2,pi/2,300), ...
    linspace(-pi/2,pi/2,300), ...
    aesa.ENUM_MODE_UNIFORM);
caxis([-30 35]);
return;
figure;
aesa.upd_plot(...
    linspace(-pi/2,pi/2,1000), ...
    0, ...
    aesa.ENUM_MODE_UNIFORM);

figure;
aesa.upd_plot(...
    0, ...
    linspace(-pi/2,pi/2,1000), ...
    aesa.ENUM_MODE_UNIFORM);