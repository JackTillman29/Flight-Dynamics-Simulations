close all;
clear all;
clc;

% Questions??
% 1 - What is the signal bandwidth?? Let's say 500kHz

%% Key parameters
Frf = mean([150 220]) * 1e6;
Fs  = 4*220e6;
kboltz = 1.38e-23;
Tref = 290;
B = Fs;
Nf = 10^(4.5/10);

% element aspect
elem_aspect = [480 384]; % from image
elem_array = [14 1];%6
elem_spacing = 1; % square spacing

arr = array('test_array');
arr.setTxProperties(Frf,1.00);
[Xe,Ye]=meshgrid(1:elem_array(1),1:elem_array(2));
arr.setElementXYZ(Xe(:),Ye(:),0*Xe(:),1);
clear Xe Ye;
arr.turnOn('all');

%% Build the interference waveform(s)
FsIF = 500e3;
DwellDistance = 300e3;
DwellTime = DwellDistance / 150 * 1e-6;
DwellSamples = floor(DwellTime*FsIF);
angOfArrival = @(dphase,lambda) asin((dphase*lambda)/(2*pi));
intfAmp = [ 1 ];
intfAz  = [ -30 18 40]*pi/180;
intfEl  = 0*[  2  15   1  6  25 60 40]*pi/180;

% Define the signal TD waveform, and where it comes from in az/el
Y = 0;
arr.setSteer(0,0);

for ksig = 1 : length(intfAmp)
    refSig = intfAmp(ksig) * single(exp(1i*2*pi*rand(1,DwellSamples)));
    refSig = refSig - mean(refSig);
    refSigAz = intfAz(ksig);
    refSigEl = intfEl(ksig);
    
    % build signals arrays
    
    
    % get array code to compute phasing to source location.
    arr.getGain(refSigAz,refSigEl);
    
    % generate elemental copies of signal
    Y = Y + (exp(1i*arr.phs).' * refSig);
end

% now apply uncorrelated random receiver noise
Y = Y + sqrt(kboltz*Tref*Nf*Fs)*single(exp(1i*2*pi*rand(elem_array(1)*elem_array(2),DwellSamples)));

% form sample-based correlation matrix
R = smpcormtx(Y);

% Form "perfect" R matrix
Rp = zeros(arr.nElements);
for k = 1 : length(intfAmp)
    arr.setSteer(intfAz(k),intfEl(k));
    vp = exp(1i*arr.phs_steer);
    Rp = Rp + ...
        vp.' * (intfAmp(k)).^2 * vp;
end
if(1)
    R = Rp;
end
[V,U]=eig(R);
totalEigVal = sum(abs(diag(U)));

figure;
plot(10*log10(abs(diag(U))),'.-');
xlabel('Eigenvalue #');
ylabel('Component Magnitude (dB)');
title('Corr. Mtx. Eig.');


%% What does it mean??
% Vectors in the columns of V seem to correspond to the spatial frequencies
% of the interference sources. Once again, they may have multiple frequency
% components in each. So **NO**, each eigenvector does not correspond to a
% single interference source.

% Steer the quiescent beam
arr.setSteer(0*pi/180,0);

% Get the quiescent steering vector
s = exp(1i*arr.phs_steer).'; % force as column vector

%arr.setSteer(0,0);

% Once this is working, an algorithm needs to be implemented that selects
% the proper number of eigenvalues to use, but for now, we know how many
% jammers there are, so we will just use that.
nEigUsed = length(intfAmp);

% Compute the inverse eigenvalue matrix (only useable terms)
Uinv = diagInv(U,nEigUsed);
Uinv(Uinv ~= 0) = 1; % each usable eigenvector will project as a unit vector


% Optimal weighting vector is R^-1 * s
w = (V * Uinv * V') * s;

arr.abf_term = (s-w)';

% here we are trying to find angles of interference
asteer = linspace(-60,60,1800);
steer_colin = 0*asteer;
azold = arr.azSteer;
elold = arr.elSteer;
for k = 1 : length(asteer)
    arr.setSteer(pi/180 * asteer(k),0);
    s = exp(1i*arr.phs_steer).';
    steer_colin(k) = abs((w'*s)/arr.nElements);
end
figure;
plot(asteer,steer_colin);
grid on;
xlabel('steer angle');
ylabel('colinearity');

% end trying to find interference angles

arr.setSteer(azold,elold);
% Compute unaltered pattern (quiescent pattern)
azval = -90:.1:90;
gOrig = arr.getAzCut(azval*pi/180,0,arr.ENUM_MODE_UNIFORM);
gAdaptive = arr.getAzCut(azval*pi/180,0,arr.ENUM_MODE_ABF);

figure;
plot(azval,20*log10(abs([gOrig' gAdaptive'])));
ylim([-50 30]);
xlim([-90 90]);
%grid on;
YL = ylim();
for k = 1 : min([length(intfAmp),nEigUsed])
    HL = line(180/pi*[intfAz(k) intfAz(k)],YL);
    set(HL,'Color','r','LineStyle',':');
end
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
title('Beamformer Response');
legend('Original','Adapted');
return;
% figure;
% FY = fft(Y,[],1);
% SFY = sum(abs(FY),2);
% plot(SFY);

figure;
plotFullPattern(arr, ...
    linspace(-pi/2,pi/2,500), ...
    linspace(-10*pi/180,pi/2,500), ...
    arr.ENUM_MODE_UNIFORM	,'azel');

hold on;
plot(180/pi*intfAz,180/pi*intfEl,'w.','MarkerSize',18);
hold off

figure;
plotFullPattern(arr, ...
    linspace(-pi/2,pi/2,500), ...
    linspace(-10*pi/180,pi/2,500), ...
    arr.ENUM_MODE_ABF	,'azel');

hold on;
plot(180/pi*intfAz,180/pi*intfEl,'w.','MarkerSize',18);