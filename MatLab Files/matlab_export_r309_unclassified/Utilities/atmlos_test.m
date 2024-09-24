close all;clear all;clc;

%% Constants and Code Setup

SPDLGT = 2.99792458e8;
M2NM =  1/1852;
REARTH = 8497000;
WAVEL = SPDLGT / 9.8e9;
xspan   = 300000;
yspan   = 200000;
xstep   = 20000;
alt_agl = 9144;
ant_agl = 5;

[X,Y,Z] = meshgrid( ...
    -xspan:xstep:xspan, ...
    -yspan:xstep:yspan, ...
    alt_agl-ant_agl);

FREQS = [.1e9,.2e9,.3e9,.6e9,1.e9,3.e9,10.e9];
ELEVS = [0.,.5,1.,2.,5.,10.];

AA = [.2739, .1881, .1605, .1031, .07371, .04119;
     .6848, .5533, .4282, .3191, .2158,  .1017;
     1.199, .9917, .7498, .5186, .3029,  .1522;
     2.210, 1.830, 1.314, .9499, .4724,  .2512;
     2.758, 2.177, 1.798, 1.168, .5732,  .3007;
     3.484, 2.592, 1.964, 1.345, .6478,  .3408;
     4.935, 3.450, 2.601, 1.718, .9130,  .4420]';

BB = [   .008648, .008644, .01106, .01723, .02313, .04076;
         .008648, .008644, .01104, .01374, .02213, .04886;
         .006837, .008795, .01110, .01474, .03116, .05360;
         .008499, .009737, .01221, .01623, .03677, .07204;
         .01030,  .01223,  .01163, .01831, .03927, .08056;
         .009745, .01225,  .01455, .02055, .04500, .08280;
         .00999,  .01340,  .01620, .02240, .03750, .08470 ]';
 


%% Round Earth computations
R  = sqrt(X.^2 + Y.^2 + Z.^2);
GR = sqrt(X.^2 + Y.^2       );

if(exist('bypass','var'))
    ELEV=ones(size(X,1),size(X,2)) .* bypass_elevation_deg;
else
    E1 = (REARTH^2+R.^2 - (REARTH+Z).^2) ./ (2*REARTH.*R);
    ELEV = 180/pi * acos(E1) - 90.0;
    ELEV(GR == 0) = 90.0;
end

RNM = min(R ./ 1852,300); % 300 nmi limit
FREQ = min(max(SPDLGT / WAVEL,FREQS(1)),FREQS(end));
ELEV = min(max(ELEV,ELEVS(1)),ELEVS(end));
FREQ_ARRAY = 0*ELEV + FREQ;

%% Interpolation Section (all input data is bounded)
% I,J indexing maps to elevation,frequency

A = interp2(FREQS,ELEVS,AA,FREQ(1:end),ELEV(1:end));
B = interp2(FREQS,ELEVS,BB,FREQ(1:end),ELEV(1:end));
A = reshape(A,size(X,1),size(X,2));
B = reshape(B,size(X,1),size(X,2));
LOSSDB = A.*(1-exp(-B.*RNM));%two way log
LOSS = 10.^(0.1*LOSSDB);%two way lin
ALOSSDB = 5*log10(LOSS);%one way log
aloss=10.^(-0.1*ALOSSDB);%one way lin
dB_per_km = ALOSSDB / R;

