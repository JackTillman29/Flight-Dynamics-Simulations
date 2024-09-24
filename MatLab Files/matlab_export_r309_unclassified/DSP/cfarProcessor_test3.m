%% 2D CFAR Example
% This example shows how to use cfarProcessor.m to process a 2D heatmap to
% detect a single target buried in noise. The kernel in this example is
% different along each dimension.

%% Define CFAR parameters
% Defining these as length-2 allows customizing each dimension in the 2D
% CFAR. Each number represents one side of the CFAR kernel in that
% dimension.

%%
% Create a CFAR kernel that has:
%%
% * dim-1,  2 guard and 3 training on one side of the cell-under-test
% * dim-2,  5 guard and 9 training on one side of the cell-under-test
nTrain = [3 9];
nGuard = [2 5];
% Given above definitions, we expect the CFAR mask size to be:
%       [2*3 + 2*2 + 1, 2*10 + 2*5 + 1] = [11,31]

%%
% Now define the CFAR threshold
cfarThreshold_dB = 13;
cfarThreshold    = 10^(cfarThreshold_dB/10);

%% Define a single cell target in noisy background
% Define the signal size and build the background noise signal matrix
N = 200;
M = 300;
x = abs(randn(N,M));
%%
% Put the target in the middle of the signal matrix
x(N/2,M/2) = 50;

%% Run the signal through the CFAR
[y,thr,cfardet] = cfarProcessor(x,cfarThreshold,nTrain,nGuard);

%% Plot the original signal, the threshold, and the detector output
figure('Position',[303 238 1136 420]);
hax(1) = subplot(1,3,1);
imagesc(x); %axis equal;

hax(2) = subplot(1,3,2);
imagesc(thr); %axis equal;

hax(3) = subplot(1,3,3);
imagesc(cfardet); %axis equal;