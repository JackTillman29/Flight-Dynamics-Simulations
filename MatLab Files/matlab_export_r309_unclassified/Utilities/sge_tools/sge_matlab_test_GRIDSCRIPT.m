% close all; clear all; clc

% xmean = 0;
% xstd  = 1;
% N = 1e5;

x = xstd*randn(1,N) + xmean;

% nbins = 100;
[xh,ed] = histcounts(x,nbins,'Normalization','pdf');

% figure; plot(xh)



