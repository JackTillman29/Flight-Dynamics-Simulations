dt = 0.1; % valid for TP
t  = 0:dt:15;
n  = length(t);

%% Build the measurement data
noiseAmplitude = 0.5;

xm = 0.0*t;
i0 = find(t < 5);
i1 = find(t >=  5.0 & t < 10);
i2 = find(t >= 10.0 & t < 15);

xm(i0) = 10.0;
xm(i1) = 0.0;
xm(i2) = 0.0 + 5.0*(t(i2)-t(i2(1)));
xtruth = xm;
xtruthd = [0 diff(xtruth)]./dt;
xtruthdd = [0 diff(xtruthd)]./dt;
xm = xm + noiseAmplitude * randn(1,n);

xHat = zeros(3,n);
xBar = zeros(3,n);
xErr = zeros(1,n);

PHI  = [ ...
    1 dt 0.5*dt*dt;
    0 1 dt;
    0 0 1 ];

%G = 1.0;
%H = 0.5;
%K = 0.2;

% TP data
G = 0.32;
H = 6.15e-2%/dt;
K = 2*2.9494e-3%/(dt*dt);

for k = 1 : n
    if(k==1)
        xBar(:,k) = [xm(k) 0 0]';
        xHat(:,k) = [xm(k) 0 0]';
        xErr(:,k) = 0.0;
    else
        % compute prediciton for this timestep
        xHat(:,k) = PHI * xBar(:,k-1);
        % compute error
        xErr(:,k) = xm(k) - xHat(1,k);
        % compute smoothed estimate
        xBar(:,k) = xHat(:,k) + [G H K]' .* xErr(:,k);
    end
end

close all;
figure(1);
subplot(3,1,1);
plot(t,xm,'ro');hold on;
plot(t,[xBar(1,:)' xtruth'],'LineWidth',2);hold off;
grid on;title('Position');legend('Measured','Filtered','Truth');
subplot(3,1,2);
plot(t,[xBar(2,:)' xtruthd'],'LineWidth',2);
grid on;title('Velocity');
subplot(3,1,3);
plot(t,[xBar(3,:)' xtruthdd'],'LineWidth',2);
grid on;xlabel('Time (sec)');
title('Acceleration');