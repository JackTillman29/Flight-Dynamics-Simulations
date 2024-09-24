function KF = StateEstimator(tdat,pdat,meas_variance,x0_std,v0_std,a0_std,Qa_std)
% This is a "quick & dirty" 9 state cartesian estimater. Consider it's use
% as "rough-in" work / a starting point.
% tdat is column time vector
% pdat is column x,y,z vector data
% This should only be used for cartesian in/out data where dt is fixed and
% only position data is known (measured)
% K. Sawmiller

%System Dynamics Matrix
A= [ ...
    0 0 0 1 0 0 0 0 0
    0 0 0 0 1 0 0 0 0
    0 0 0 0 0 1 0 0 0
    0 0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 1 0
    0 0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 ];

H = [ ...
    1 0 0 0 0 0 0 0 0
    0 1 0 0 0 0 0 0 0
    0 0 1 0 0 0 0 0 0];

%Measurement covariance
Rk = [ ...
    1 0 0
    0 1 0
    0 0 1 ] * meas_variance;

%Initial state covariance
Pk = diag([x0_std^2 x0_std^2 x0_std^2 v0_std^2 v0_std^2 v0_std^2 a0_std.^2 a0_std.^2 a0_std.^2 ]);
I= eye(9);
% Process Noise
Q = diag([0 0 0 0 0 0 1 1 1]*Qa_std.^2);
% compute discrete process noise and state transition matrix
dt_mean = mean(diff(tdat));

Phi = expm(A.*dt_mean);
Qk = ComputeNumericQk(A,Q,dt_mean,10);


N = length(tdat);

KF.xhat = zeros(N,9);
KF.xcov = zeros(N,9);

xhat = [pdat(1,:) (pdat(2,:)-pdat(1,:))/dt_mean 0 0 0]';
KF.xhat(1,:) = xhat;
KF.err = zeros(N,3);

for krow = 2 : N
    %Solve Riccati Equation for Mk (pre-measurement covariance)
    Mk= Phi * Pk * Phi' + Qk;
    
    % Solve Riccati Equation for Kk (Kalman Gains)
    Kk = Mk * H' * inv(H * Mk * H' + Rk);
    
    %Solve Riccati Equation for Pk (post-measurement covariance)
    Pk =(I - Kk * H) * Mk;
    
    %Compute expected measurement at this timestep (due to non-linear H matrix)
    xbar = Phi * xhat; %last estimate propagated to now
    
    % Compute expected measurement
    xm_expected = H * xbar;
    % Get current measurement
    xm = pdat(krow,:)';
    % Compute residual for this timestep (measured -expected)
    residual= xm - xm_expected;
    KF.err(krow,:) = residual';
    
    % Compute new state estimate
    xhat = xbar + Kk * residual ;
    %Save data
    KF.xhat(krow,:) = xhat';
    KF.xcov(krow,:) = Pk(1:(9+1):(9^2));
end
