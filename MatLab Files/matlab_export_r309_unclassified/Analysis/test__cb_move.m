close all; clear all; clc


figure('Position',[205 234 560 420]);
for k= 1:4
subplot(2,2,k)
plot(rand*randn(1,100))
end

figure('Position',[771 209 560 420]);
add_analysis_callbacks;
plot(randn(1,50));






