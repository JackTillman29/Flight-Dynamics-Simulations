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
elem_array = [14 6];%6
elem_spacing = 1; % square spacing

arr = array('test_array');
arr.setTxProperties(Frf,1.00);
[Xe,Ye]=meshgrid(1:elem_array(1),1:elem_array(2));
arr.setElementXYZ(Xe(:),Ye(:),0*Xe(:),1);
clear Xe Ye;
arr.turnOn('all');


if(0)
    g=arr.getGainTable( ...
        -pi/180 * [-45:45], ...
        -pi/180 * [-30:30], ...
        arr.ENUM_MODE_UNIFORM);
end
%% Build the interference waveform(s)
FsIF = 500e3;
DwellDistance = 300e3;
DwellTime = DwellDistance / 150 * 1e-6;
DwellSamples = floor(DwellTime*FsIF);
angOfArrival = @(dphase,lambda) asin((dphase*lambda)/(2*pi));
intfAmp = [ 1  1 1 1];
intfAz  = [ -30 -10 20 40]*pi/180;
intfEl  = [ 0 0 0 6]*pi/180;

% Define the signal TD waveform, and where it comes from in az/el
Y = 0;
for ksig = 1 : length(intfAmp)
    refSig = intfAmp(ksig) * single(exp(1i*2*pi*rand(1,DwellSamples)));
    refSig = refSig - mean(refSig);
    refSigAz = intfAz(ksig);
    refSigEl = intfEl(ksig);
    
    % build signals arrays
    
    arr.setSteer(0,0);
    % get array code to compute phasing to source location.
    arr.getGain(refSigAz,refSigEl);
    
    % generate elemental copies of signal
    Y = Y + (exp(1i*arr.phs).' * refSig);
end

% now apply uncorrelated random receiver noise
Y = Y + sqrt(kboltz*Tref*Nf*Fs)*single(exp(1i*2*pi*rand(elem_array(1)*elem_array(2),DwellSamples)));

% form sample-based correlation matrix
R = smpcormtx(Y);

[V,U]=eig(R);
totalEigVal = sum(abs(diag(U)));

figure;
plot(10*log10(abs(diag(U))),'.-');
xlabel('Eigenvalue #');
ylabel('Component Magnitude (dB)');
title('Corr. Mtx. Eig.');
% eigValThresh = 1e-5;
% eigAz = 0*U(:,1);
% eigEl = 0*U(:,1);
% for kv = 1 : size(U,1)
%     % if eigenvalue is deemed significant
%     if(abs(U(kv,kv)) > eigValThresh)
%         % pick off eigenvector, reshape into "array" space (2D)
%         thisEV = reshape(V(:,kv),elem_array(2),elem_array(1));
%         % compute horizontal accumulated phase
%         accum_phase = unwrap(angle(thisEV(1,:)));
%         phase_per_element = (accum_phase(end) - accum_phase(1))/(elem_array(1)-1);
%         eigAz(kv) = angOfArrival(phase_per_element,arr.wavelength)*180/pi;
%         text(kv,double(10*log10(abs(U(kv,kv)))),num2str(eigAz(kv)),'HorizontalAlignment','center');
%     else
%         % do nothing.
%         %disp(sprintf('Eigenvalue %8.3e skipped',abs(U(kv,kv))));
%     end
% end

%% What does it mean??
% Vectors in the columns of V seem to correspond to the spatial frequencies
% of the interference sources. Once again, they may have multiple frequency
% components in each. So **NO**, each eigenvector does not correspond to a
% single interference source.

% Steer the quiescent beam
arr.setSteer(0,0);

% Get the quiescent steering vector
s = exp(1i*2*pi*arr.phs_steer);

% Once this is working, an algorithm needs to be implemented that selects
% the proper number of eigenvalues to use, but for now, we know how many
% jammers there are, so we will just use that.
nEigUsed = length(intfAmp);

% Compute the inverse eigenvalue matrix (only useable terms)
Uinv = diagInv(U,nEigUsed);
Uinv(Uinv ~= 0) = 1; % each usable eigenvector will project as a unit vector


% Optimal weighting vector is R^-1 * s
w = (V * Uinv * V') * s.';

%w = computeAdaptiveWeight(s,diag(U),V,'direct',nEigUsed); % WORKS


arr.abf_term = s-w';

%eig_pinv = Vr * UrInv * Vr'; % only because eig(A * A'), so eigvecs are orthonormal
%arr.abf_term =  s-((Vr(:,1)'*s')*Vr(:,1))'; -- works, but no eigenvalue???
%arr.abf_term = s-(Vr*Vr' * s')';

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

%return
figure;
plotFullPattern(arr,linspace(-pi,pi,500),linspace(-pi/2,pi/2,250),arr.ENUM_MODE_ABF	,'azel');