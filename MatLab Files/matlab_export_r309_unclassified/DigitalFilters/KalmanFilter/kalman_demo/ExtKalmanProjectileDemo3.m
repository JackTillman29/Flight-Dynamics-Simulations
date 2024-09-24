% Script to show how a Kalman filter works and is benificial
clear;clc;
% First, generate the data set for the analysis.  The script will generate
% both the truth data and measured data for the duration of the flight. It
% calculates the impact time and only computes data to that point.

% Specify Time Step
dt = 0.1;

% Specify sensor properties
radar_measurement_error_range = 50;    % meters
radar_measurement_error_angle = 3*0.0175; % radians
radar_position_x = 85;
radar_position_y = 0;

% Specify Cannon Initial Conditions
CannonPosition_x = 0;
CannonPosition_y = 0;
MuzzleVelocity_x = 2500;
MuzzleVelocity_y = 2250;

% Generate Data
% [t,TruthData,MeasuredData]=GenerateData( ...
%     CannonPosition_x,CannonPosition_y, ...
%     MuzzleVelocity_x,MuzzleVelocity_y, ...
%     radar_position_x,radar_position_y, ...
%     radar_measurement_error_range,radar_measurement_error_angle, ...
%     dt);

[t,TruthData,MeasuredData]=GenerateData2( ...
    CannonPosition_x,CannonPosition_y, ...
    MuzzleVelocity_x,MuzzleVelocity_y, ...
    radar_position_x,radar_position_y, ...
    radar_measurement_error_range,radar_measurement_error_angle, ...
    dt,1e-4);


% Define the non-linear H matrix function
ComputeH =@(xt,yt,xs,ys) ...
    [ ...
    (xt-xs)/sqrt((yt-ys)^2+(xt-xs)^2) ...
    0 ...
    (yt-ys)/sqrt((yt-ys)^2+(xt-xs)^2) ...
    0;
    -(yt-ys)/((yt-ys)^2+(xt-xs)^2) ...
    0 ...
    (xt-xs)/((yt-ys)^2+(xt-xs)^2) ...
    0 ...
    ];

% Define the initial state estimate
xhat = [ ...
    0
    100
    0
    100
    ];

QX = 10^3;
QY = 10^3;

% Define the initial state estimate covariance
Pk = [ ...
    50^2       0       0       0
    0       1e12       0       0
    0       0       50^2       0
    0       0       0       1e12
    ];

% Define the sensor measurement covariance
Rk = [ ...
    radar_measurement_error_range^2     0
    0                                   radar_measurement_error_angle^2
    ];

% Define the state transition matrix, Phi
Phi = [
    1  dt   0   0
    0   1   0   0
    0   0   1  dt
    0   0   0   1
    ];

% Define the discrete process noise matrix (integrate(phi*Q*phi',dt))
Qk = [ ...
    QX*dt^3/3   QX*dt^2/2       0           0
    QX*dt^2/2   QX*dt           0           0
    0           0               QY*dt^3/3   QY*dt^2/2
    0           0               QY*dt^2/2   QY*dt ...
    ];

% Define the feedforward term due to gravity (integrate(Phi*G,dt))
Gk = [
    0
    0
    -9.80665*dt^2/2
    -9.80665*dt
    ];

% Define the Identity Matrix for this filter
I = eye(4);

% Create Data Storage for Post Processing
KF.xhat = zeros(length(t),4);   % state estimates
KF.xcov = KF.xhat;              % state covariance estimates
KF.xerr = KF.xhat;              % truth-state errors


% Run Filter
KF.xhat(1,:) = xhat';
KF.xcov(1,:) = Pk([1 6 11 16]);


for k = 2:length(t)
    
    % Define measurement matrix for this timestep
    xm = [MeasuredData.range(k) MeasuredData.angle(k)]';
    
    % Define the H matrix for this timestep
    H = ComputeH(xhat(1),xhat(3),radar_position_x,radar_position_y);

    % Solve Riccati Equation for Mk (pre-measurement covariance)
    Mk = Phi * Pk * Phi' + Qk;

    % Solve Riccati Equation for Kk (Kalman Gains)
    Kk = Mk * H' * inv(H * Mk * H' + Rk);

    % Solve Riccati Equation for Pk (post-measurement covariance)
    Pk = (I - Kk * H) * Mk;

    % Compute expected measurement at this timestep (due to non-linear H
    % matrix)
    xbar = Phi * xhat + Gk;  %last estimate propagated to now + gravity term
    xm_expected = [ ...
        sqrt( (xbar(1)-radar_position_x)^2 + (xbar(3)-radar_position_y)^2)
        atan2((xbar(3)-radar_position_y),(xbar(1)-radar_position_x))];
    
    % Compute residual for this timestep (measured - expected)
    residual = xm - xm_expected;
    
    % Compute new state estimate
    xhat = xbar + Kk * residual ;
    
    % Save data
    KF.xhat(k,:) = xhat';
    KF.xcov(k,:) = Pk([1 6 11 16]);
end

KF.xerr =KF.xhat - [TruthData.position_x TruthData.velocity_x TruthData.position_y TruthData.velocity_y];


MeasuredX = MeasuredData.range .* cos(MeasuredData.angle);
MeasuredY = MeasuredData.range .* sin(MeasuredData.angle);


figure(1);
plot( MeasuredX, MeasuredY, '.', ...
    TruthData.position_x, TruthData.position_y,'.');
grid on;
xlabel('Downrange Meters');
ylabel('Crossrange Meters');
legend('Measured Data','Truth Data');

figure(2);
subplot(2,1,1);
plot(t(2:end),[diff(MeasuredX)./dt diff(TruthData.position_x)./dt]);
grid on;
xlabel('Time (sec)');
ylabel('X Velocity (m/s)');
legend('Dead Reckoned','True Velocity');
subplot(2,1,2);
plot(t(2:end),[diff(MeasuredY)./dt diff(TruthData.position_y)./dt]);
grid on;
xlabel('Time (sec)');
ylabel('Y Velocity (m/s)');

figure(3);
plot( KF.xhat(:,1), KF.xhat(:,3), '.', ...
    TruthData.position_x, TruthData.position_y,'.');
grid on;
xlabel('Downrange Meters');
ylabel('Crossrange Meters');
legend('Kalman Estimates','Truth Data');

figure(4);
subplot(2,1,1);
plot(t,[KF.xhat(:,2) TruthData.velocity_x]);
grid on;
xlabel('Time (sec)');
ylabel('X Velocity (m/s)');
legend('Kalman Estimate','True Velocity');
subplot(2,1,2);
plot(t,[KF.xhat(:,4) TruthData.velocity_y]);
grid on;
xlabel('Time (sec)');
ylabel('Y Velocity (m/s)');

figure(5);
subplot(2,1,1);
plot( ...
    t,sqrt(KF.xcov(:,1)),'r', ...
    t,-sqrt(KF.xcov(:,1)),'r', ...
    t,KF.xerr(:,1),'b');
grid on;
xlabel('Time (sec)');
ylabel('X Deviation');
subplot(2,1,2);
plot( ...
    t,sqrt(KF.xcov(:,3)),'r', ...
    t,-sqrt(KF.xcov(:,3)),'r', ...
    t,KF.xerr(:,3),'b');
grid on;
xlabel('Time (sec)');
ylabel('Y Deviation');

figure(6);
subplot(2,1,1);
plot( ...
    t,sqrt(KF.xcov(:,2)),'r', ...
    t,-sqrt(KF.xcov(:,2)),'r', ...
    t,KF.xerr(:,2),'b');
grid on;
xlabel('Time (sec)');
ylabel('X Vel Deviation');
subplot(2,1,2);
plot( ...
    t,sqrt(KF.xcov(:,4)),'r', ...
    t,-sqrt(KF.xcov(:,4)),'r', ...
    t,KF.xerr(:,4),'b');
grid on;
xlabel('Time (sec)');
ylabel('Y Vel Deviation');

%close([1 2 4 5 6]);