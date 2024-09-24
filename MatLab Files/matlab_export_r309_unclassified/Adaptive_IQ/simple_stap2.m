close all;
clear all;
clc;

addpath('Z:\sawmillkd\MATLAB/DSP');
addpath('Z:\sawmillkd\MATLAB/Windows');
addpath('Z:\sawmillkd\MATLAB/Printing');
%pp = PrepForPrint;

% basic array properties
N       = 11;
fTx     = 10e9;
c       = 3.0e8;
lambda  = c / fTx;
d2r     = pi / 180;

% SNR

% element spacing
d = 0.5 * lambda;

% Taper
w_taper = taylorwin(N,5,-50);

% steering vector
theta0 = 0*d2r;

% plot properties
sinspace =0;

calc_s = inline('exp(1i*2*pi*((0:(N-1))'')*d/lambda*sin(theta0))');
calc_pwr_pat = inline('20*log10(abs(fftshift(fft(weights.*s,NFFT))))','weights','s','NFFT');
wgt_norm = inline('x./sqrt(x''*x)');

s = calc_s(N,.5*lambda,lambda,theta0);

NFFT = 8192;
nrmang = linspace(-d/lambda,d/lambda,NFFT);
angvec =  asin(lambda/d * nrmang) ./ d2r;

figure;
if(sinspace)
    xvec = nrmang;
    xlab = '(d/\lambda) sin(\theta)';
else
    xvec = angvec;
    xlab = '\theta (deg)';
end
plot(xvec,calc_pwr_pat(wgt_norm(w_taper),s,NFFT));
%ylim([-50 0]);
xlabel(xlab);ylabel('Relative Power (dB)');grid on;
title('Beamformer Response');
add_print_callbacks;



%% Case 1: Additive White Noise
const_k = 1.38e-23;
Nf      = 10^(5/10);
Tref    = 290;
B       = 160e3;
sigma2  = Nf * const_k * Tref * B;

Rthermal = diag(ones(1,N)*sigma2);

%% Case 2: Colored Noise w/ Jammer
% ===============================================
%
% NOTE: THIS SECTION PLACES POWER IN MULTIPLES OF
%       NOISE FLOOR.
%
% ===============================================
angJams = [-40 40]*d2r;
sigmaJs = [10 10];
R = 0;
Rthermal = eye(N) * 1;
for k = 1 : length(angJams)
    sj = calc_s(N,.5*lambda,lambda,angJams(k));
    sigmaJ = sigmaJs(k);
    R = R + sigmaJ * (sj * sj');
end
R = R + Rthermal;

% Generate eigenspace data R = v * u * v'
[v,u] = eig( R );
eigvals = diag(u);
eigvals_abs = abs(eigvals);


%min_eig = sigma2;
min_eig = 2;

% compute weighted sum "principal components" method
wsum = 0;
s_tgt = calc_s(N,0.5*lambda,lambda,theta0);

for k = 1 : N
    % pull out unit vector associated w/ this eigenvalue
    un = v(:,k);
    
    wsum = wsum + ...
        (max(eigvals_abs(k)-min_eig,0))/eigvals_abs(k) * (un' * s_tgt ) * un;
end
w_eig = s_tgt-wsum;



figure;
plot(xvec,[ ...
    calc_pwr_pat(wgt_norm(w_taper),s,NFFT) ...
    calc_pwr_pat(wgt_norm(w_taper.*w_eig),s,NFFT) ...
    ]);
add_print_callbacks;
xlabel(xlab);
ylabel('Gain (dBi)');
title('Adapted Beamformer');grid on;
legend('Normal','PrincipalComp'); 


%% Generate Real Data
tgtAz  = 30*d2r;
tgtPwr = 0; % power at receiver
jamAz  = -7.6*d2r;
jamPwr = 1*1000; % power at receiver (factor x noise floor)

% Populate with receiver noise
%s = sqrt(sigma2) * exp(1i*2*pi*randn(N,1));

% Add target signal
w_tgt = calc_s(N,.5*lambda,lambda,tgtAz)./(1/sqrt(N));
%s = s + sqrt(tgtPwr) * calc_s(N,.5*lambda,lambda,tgtAz);

% Add jamming signal
%s = s + sqrt(jamPwr) * calc_s(N,.5*lambda,lambda,jamAz);

% Get Target Return -- TODO: Check array math
%signal_target = abs(w_tgt' * s);

%
figure;
%subplot(2,1,1);
plot(xvec,calc_pwr_pat(wgt_norm(w_taper),w_tgt,NFFT));hold on;
plot(tgtAz/d2r,20,'bv','MarkerFaceColor','b');
plot(jamAz/d2r,20,'rv','MarkerFaceColor','r');
xlabel(xlab);
ylabel('Power (dBW)');
title('Steered Non-Adapted Beamformer');grid on;
ylim([-50 30]);

%% Synthetically generate sample data (no target)
L = 10; % number of times you sample the elements
R_SMI = 0;

for k = 1 : L
    % noise
    jam_rand_phase  = (exp(1i*2*pi*rand(1)) * ones(1,N)).';
    jam_rand_phase2 = (exp(1i*2*pi*rand(1)) * ones(1,N)).';
    
    jam_dir_phase  = calc_s(N,.5*lambda,lambda,jamAz);
    s_smi = 1 * jam_rand_phase + sqrt(jamPwr) .* jam_rand_phase2 .* jam_dir_phase;
    
    R_SMI = R_SMI + s_smi*s_smi';
end
R_SMI = R_SMI ./ L;

%w_jam = calc_s(N,.5*lambda,lambda,jamAz);
%R = jamPwr * (w_jam * w_jam') + eye(N);
%w_pfc = R\(w_tgt.*w_taper); % inv(R) * (sw)

% now apply principal component analysis
% Generate eigenspace data R = v * u * v'
[v,u] = eig( R_SMI );
eigvals = diag(u);
eigvals_abs = abs(eigvals);
eigvals_pca = 0*eigvals_abs;
min_eig = 2;

% compute weighted sum "principal components" method
wsum = 0;

for k = 1 : N
    % pull out unit vector associated w/ this eigenvalue
    un = v(:,k);
    eigvals_pca(k) = (max(eigvals_abs(k)-min_eig,0))/eigvals_abs(k);
    wsum = wsum + ...
        eigvals_pca(k) * (un' * w_tgt ) * un;
end
w_eig = w_tgt-wsum;



% Update Beamformer
%w_smi = R_SMI\(w_tgt.*w_taper); % inv(R) * (sw)

figure;
plot(xvec,[ ...
    calc_pwr_pat(wgt_norm(w_taper.*w_tgt),0*w_tgt+1,NFFT), ...
    calc_pwr_pat(wgt_norm(w_taper.*w_eig),0*w_tgt+1,NFFT) ...
    calc_pwr_pat(w_eig,0*w_tgt+1,NFFT) ...
    ]);hold on;
plot(tgtAz/d2r,20,'bv','MarkerFaceColor','b');
plot(jamAz/d2r,20,'rv','MarkerFaceColor','r');
xlabel(xlab);
ylabel('Power (dBW)');
title('Adapted Beamformer');grid on;
ylim([-50 30]);
legend('Non-Adaptive','PCA');

return;
%% Jammer Loop Animation
vid = 0;
if(vid==1)
    vidobj = VideoWriter('jamloop.avi');
    vidobj.Quality = 95;
    vidobj.FrameRate = 30;
end
jamAzList = -90:1:90;
temp=(0:(length(jamAzList)-1));
temp = temp ./ temp(end);
tgtAzStart = 10;
tgtAzEnd   = 10;
tgtAzDif   = tgtAzEnd - tgtAzStart;

for kjam = 1 : length(jamAzList)

%tgtAz  = 30*d2r;
tgtAz   = (tgtAzStart + temp(kjam)*tgtAzDif)*d2r;
tgtPwr = 0; % power at receiver
jamAz  = jamAzList(kjam)*d2r;
jamPwr = 50000; % power at receiver (factor x noise floor)

% Populate with receiver noise
s = sqrt(sigma2) * exp(1i*2*pi*randn(N,1));

% Add target signal
w_tgt = calc_s(N,.5*lambda,lambda,tgtAz)./(1/sqrt(N));
s = s + sqrt(tgtPwr) * calc_s(N,.5*lambda,lambda,tgtAz);

% Add jamming signal
s = s + sqrt(jamPwr) * calc_s(N,.5*lambda,lambda,jamAz);

% Get Target Return -- TODO: Check array math
signal_target = abs(w_tgt' * s);

% Synthetically generate sample data (no target)
L = 2*N;
R_SMI = 0;
for k = 1 : L
    % noise
    jam_rand_phase  = (exp(1i*2*pi*rand(1)) * ones(1,N)).';
    jam_rand_phase2 = (exp(1i*2*pi*rand(1)) * ones(1,N)).';
    jam_dir_phase  = calc_s(N,.5*lambda,lambda,jamAz);
    
    s_smi = 1 * jam_rand_phase + sqrt(jamPwr) .* jam_rand_phase2 .* jam_dir_phase;
    R_SMI = R_SMI + s_smi*s_smi';
end
R_SMI = R_SMI ./ L;

% now apply principal component analysis
% Generate eigenspace data R = v * u * v'
[v,u] = eig( R_SMI );
eigvals = diag(u);
eigvals_abs = abs(eigvals);

min_eig = 2;

% compute weighted sum "principal components" method
wsum = 0;

for k = 1 : N
    % pull out unit vector associated w/ this eigenvalue
    un = v(:,k);
    
    wsum = wsum + ...
        (max(eigvals_abs(k)-min_eig,0))/eigvals_abs(k) * (un' * w_tgt ) * un;
end
w_eig = w_tgt-wsum;


if(kjam == 1)
    figure;
    hp=plot(xvec, ...
        calc_pwr_pat(wgt_norm(w_taper.*w_eig),0*w_tgt+1,NFFT) ...
        );
    hold on;
    ht=plot(tgtAz/d2r,20,'bv','MarkerFaceColor','b');
    hj=plot(jamAz/d2r,20,'rv','MarkerFaceColor','r');
    xlabel(xlab);
    ylabel('Antenna Gain (dBi)');
    title('Adapted Beamformer');grid on;
    ylim([-50 30]);
    if(vid)
        open(vidobj);
        disp('set up screen then press any key');
        add_print_callbacks;
        pause;
    end
    
else
    set(hp,'YData',calc_pwr_pat(wgt_norm(w_taper.*w_eig),0*w_tgt+1,NFFT));
    set(hj,'XData',jamAz/d2r);
    set(ht,'XData',tgtAz/d2r);
end
drawnow;
if(vid)
     writeVideo(vidobj,getframe(gcf));
end
end
close(vidobj);
%% Multiple Jammers Animation
npoints = 900;

vid = 0;
if(vid==1)
    vidobj = VideoWriter('multijamloop.avi');
    vidobj.Quality = 95;
    vidobj.FrameRate = 30;
end

temp=(0:(npoints-1));
temp = temp ./ temp(end);
tgtAzStart = 0;
tgtAzEnd   = 0;
tgtAzDif   = tgtAzEnd - tgtAzStart;
jamAzStart = [-60 -50 10 20];
jamAzEnd   = [-10 -40 30 80];
jamAzDif   = jamAzEnd - jamAzStart;
clear jamAz;
for kjam = 1 : npoints

%tgtAz  = 30*d2r;
tgtAz   = (tgtAzStart + temp(kjam)*tgtAzDif)*d2r;
tgtPwr = 0; % power at receiver
jamAzVec  = (jamAzStart + temp(kjam)*jamAzDif)*d2r;
jamPwr = 50000; % power at receiver (factor x noise floor)

% Populate with receiver noise
s = sqrt(sigma2) * exp(1i*2*pi*randn(N,1));

% Add target signal
w_tgt = calc_s(N,.5*lambda,lambda,tgtAz)./(1/sqrt(N));
s = s + sqrt(tgtPwr) * calc_s(N,.5*lambda,lambda,tgtAz);

% Get Target Return -- TODO: Check array math
signal_target = abs(w_tgt' * s);

% Synthetically generate sample data (no target)
L = 2*N;
R_SMI = 0;
for k = 1 : L
    % noise
    noise_rand_phase  = (exp(1i*2*pi*rand(1)) * ones(1,N)).';
    
    
    
    
     s_smi = 1 * noise_rand_phase;
    
    for ijammer = 1 : length(jamAzStart)
        jam_rand_phase = (exp(1i*2*pi*rand(1)) * ones(1,N)).';
        jam_dir_phase  = calc_s(N,.5*lambda,lambda,jamAzVec(ijammer));
        s_smi = s_smi + sqrt(jamPwr) .* jam_rand_phase .* jam_dir_phase;
    end
    R_SMI = R_SMI + s_smi*s_smi';
end
R_SMI = R_SMI ./ L;

% now apply principal component analysis
% Generate eigenspace data R = v * u * v'
[v,u] = eig( R_SMI );
eigvals = diag(u);
eigvals_abs = abs(eigvals);

min_eig = 2;

% compute weighted sum "principal components" method
wsum = 0;

for k = 1 : N
    % pull out unit vector associated w/ this eigenvalue
    un = v(:,k);
    
    wsum = wsum + ...
        (max(eigvals_abs(k)-min_eig,0))/eigvals_abs(k) * (un' * w_tgt ) * un;
end
w_eig = w_tgt-wsum;


if(kjam == 1)
    figure;
    hp=plot(xvec, ...
        calc_pwr_pat(wgt_norm(w_taper.*w_eig),0*w_tgt+1,NFFT) ...
        );
    hold on;
    ht=plot(tgtAz/d2r,20,'bv','MarkerFaceColor','b');
    hj = [];
    for ijammer = 1 : length(jamAzStart)
        hj=[hj plot(jamAzVec(ijammer)./d2r,20,'rv','MarkerFaceColor','r')];
    end
    
    xlabel(xlab);
    ylabel('Antenna Gain (dBi)');
    title('Adapted Beamformer');grid on;
    ylim([-50 30]);
    xlim([-90 90]);
    if(vid)
        open(vidobj);
        disp('set up screen then press any key');
        add_print_callbacks;
        pause;
    end
else
    set(hp,'YData',calc_pwr_pat(wgt_norm(w_taper.*w_eig),0*w_tgt+1,NFFT));
    for ijammer = 1 : length(jamAzStart)
       set(hj(ijammer),'XData',jamAzVec(ijammer)./d2r);
    end
    
    set(ht,'XData',tgtAz/d2r);
end
drawnow;
if(vid)
     writeVideo(vidobj,getframe(gcf));
end
end