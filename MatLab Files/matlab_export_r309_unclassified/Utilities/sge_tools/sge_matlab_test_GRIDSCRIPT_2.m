

% 1-watt, additive white guassian noise generator
awgn          = @(N,M) sqrt(0.5)*( randn(N,M) + 1i*randn(N,M) );
calcMeanPwr_W = @(x,z0) mean(abs(x).^2/z0,2);
w2dbm         = @(x) 10*log10(x) + 30;

% Rout = 50;
% 
% N     = 1e6;
% nreps = 5;

% knob_avgPwr_1_W = 1e-3;
% knob_avgPwr_2_W = 1e-3;

% knob_avgPwr_1_dBm_arr = [-30:1:0];
% knob_avgPwr_1_dBm_arr = [-30:1:0];

% n1 = length(knob_avgPwr_1_dBm_arr);
% n2 = length(knob_avgPwr_2_dBm_arr);
% meanPwr_dBm = zeros(n1,n2)

% k1 = 0;
% for knob_avgPwr_1_dBm = knob_avgPwr_1_dBm_arr
%     k1 = k1 + 1;
%     k2 = 0;
% for knob_avgPwr_2_dBm = knob_avgPwr_1_dBm_arr
%     k2 = k2 + 1;
    
knob_avgPwr_1_W = 10^((knob_avgPwr_1_dBm-30)/10);
knob_avgPwr_2_W = 10^((knob_avgPwr_2_dBm-30)/10);

x1 = sqrt(Rout*knob_avgPwr_1_W)*awgn(nreps,N);
x2 = sqrt(Rout*knob_avgPwr_2_W)*awgn(nreps,N);

x = x1 + x2;

meanPwr_W = calcMeanPwr_W(x,Rout);

% compute average of all reps
meanPwr_W = mean(meanPwr_W);

meanPwr_dBm(k1,k2) = w2dbm(meanPwr_W);

% clean up workspace so that a ton of data isn't saved into matfiles
clear meanPwr_W x x1 x2

% end
% end







