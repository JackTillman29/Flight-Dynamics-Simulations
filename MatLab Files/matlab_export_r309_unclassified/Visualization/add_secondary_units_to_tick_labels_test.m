close all; clear all; clc

df = 123*1e-3;
f_offset = 0;

N = 100;
x = -ceil(N/2):(floor(N/2)-1);
% x = x .* df
y = randn(1,N);

figure;
add_print_callbacks;
plot(x,y);

xlabel('Doppler Bin')

% add_secondary_units_to_tick_labels(gca,'x',@(x) round(x/df),'bin')
add_secondary_units_to_tick_labels(gca,'x',@(x) x*df,'kHz')


set(gca,'FontSize',14)
set(gca,'FontWeight','bold')
set(gca,'XTickLabelRotation',45)



% msgbox(['Now zoom in, pan around the window, change xticks, etc... ' ...
%     'then right-click in the axes and select ''update secondary tick labels'''])






