close all;
clear all;
clc;

%addpath('Z:\sawmillkd\MATLAB/DSP');
%addpath('Z:\sawmillkd\MATLAB/Windows');
%addpath('Z:\sawmillkd\MATLAB/Printing');
%pp = PrepForPrint;

% basic array properties
N       = 14;
fTx     = 185e6;
c       = 3.0e8;
lambda  = c / fTx;
d2r     = pi / 180;

% SNR

% element spacing
d = 1;

% Taper
w_taper = taylorwin(N,5,-50);
w_taper = 0*w_taper + 1;

% steering vector
theta0 = 20*d2r;

% plot properties
sinspace =1;

calc_s = inline('exp(1i*2*pi*((0:(N-1))'')*d/lambda*sin(theta0))');
calc_pwr_pat = inline('20*log10(abs(fftshift(fft(weights.*s,NFFT))))','weights','s','NFFT');
wgt_norm = inline('x./sqrt(x''*x)');

s = calc_s(N,d,lambda,theta0);

NFFT = 8192;
nrmang = linspace(-1,1,NFFT);
angvec =  asin(nrmang) ./ d2r;

figure;
if(sinspace)
    xvec = nrmang;
    xlab = '(d/\lambda) sin(\theta)';
else
    xvec = angvec;
    xlab = '\theta (deg)';
end
%plot(xvec,calc_pwr_pat(wgt_norm(w_taper),s,NFFT));
[p,a]=calc_pwr_pat2(wgt_norm(w_taper).*s,d,lambda);
%ylim([-50 0]);
plot(a,p);
xlabel('Deg');ylabel('Relative Power (dB)');grid on;
title('Beamformer Response');
add_print_callbacks;
ylim([-80 20]);

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
angJams = [-10 40]*d2r;
sigmaJs = [10 10];
R = 0;
Rthermal = eye(N) * 1;
for k = 1 : length(angJams)
    sj = calc_s(N,d,lambda,angJams(k));
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
s_tgt = calc_s(N,d,lambda,theta0);


for k = 1 : N
    % pull out unit vector associated w/ this eigenvalue
    un = v(:,k);
    
    wsum = wsum + ...
        (max(eigvals_abs(k)-min_eig,0))/eigvals_abs(k) * (un' * s_tgt ) * un;
end
w_eig = s_tgt-wsum;


% ============ attempt keith's direct approach
if(1)
    Uinv = diagInv(u,length(angJams));
    Uinv(Uinv ~= 0) = 1; % each usable eigenvector will project as a unit vector

    w = (v * Uinv * v') * s_tgt;

    %w = computeAdaptiveWeight(conj(s),diag(u),v,'direct',length(angJams))'; % WORKS
    w_eig = (s_tgt-w);
    
end


[p_orig,a]=calc_pwr_pat2(wgt_norm(w_taper).*s,d,lambda);
[p_abf]=calc_pwr_pat2(wgt_norm(w_taper).*w_eig,d,lambda);
figure;
plot(a,[ ...
    p_orig.' ...
    p_abf.' ...
    ]);
add_print_callbacks;
xlabel('degrees');
ylabel('Gain (dBi)');
title('Adapted Beamformer');grid on;
legend('Normal','PrincipalComp'); 
ylim([-80 20]);


