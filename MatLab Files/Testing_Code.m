%Booz Allen Hamilton
%Jack Tillman

%MatLab testing script - this script is used to test code and syntax

% Measure runtime
tic

%%Simple Harmonic Motion Analysis

% Parameters
m = 5;       % Mass
k = 10;       % Spring constant
omega = sqrt(k/m); % Angular frequency

% Time settings
t0 = 0;      % Start time
tf = 25;     % End time
dt = 0.1; % Reduced time step to quadruple the number of time points
N = (tf - t0) / dt; % Number of steps

% Initial conditions
x0 = 1;      % Initial displacement
v0 = 0;      % Initial velocity
y0 = [x0; v0]; % Initial state vector

% Define the time vector
t = t0:dt:tf;

% Analytical solution
A = x0; % Amplitude
analytical_sol = A * cos(omega * t);

% Preallocate result arrays
euler_sol = zeros(2, length(t));
heun_sol = zeros(2, length(t));
rk4_sol = zeros(2, length(t));

euler_sol(:,1) = y0;
heun_sol(:,1) = y0;
rk4_sol(:,1) = y0;

% Euler Method Integration
for i = 1:length(t)-1
    euler_sol(:,i+1) = euler_sol(:,i) + dt * f(t(i), euler_sol(:,i), m, k);
end

% Heun's Method Integration
for i = 1:length(t)-1
    y = heun_sol(:,i);
    k1 = dt * f(t(i), y, m, k);
    k2 = dt * f(t(i) + dt, y + k1, m, k);
    heun_sol(:,i+1) = y + (k1 + k2) / 2;
end

% RK4 Integration
for i = 1:length(t)-1
    k1 = dt * f(t(i), rk4_sol(:,i), m, k);
    k2 = dt * f(t(i) + dt/2, rk4_sol(:,i) + k1/2, m, k);
    k3 = dt * f(t(i) + dt/2, rk4_sol(:,i) + k2/2, m, k);
    k4 = dt * f(t(i) + dt, rk4_sol(:,i) + k3, m, k);
    
    rk4_sol(:,i+1) = rk4_sol(:,i) + (k1 + 2*k2 + 2*k3 + k4) / 6;
end

% Compute differences between numerical solutions and analytical solution
euler_diff = abs(euler_sol(1,:) - analytical_sol);
heun_diff = abs(heun_sol(1,:) - analytical_sol);
rk4_diff = abs(rk4_sol(1,:) - analytical_sol);

% Differences between numerical integrators
euler_vs_heun = abs(euler_sol(1,:) - heun_sol(1,:));
euler_vs_rk4 = abs(euler_sol(1,:) - rk4_sol(1,:));
heun_vs_rk4 = abs(heun_sol(1,:) - rk4_sol(1,:));

% Plot results
figure;

% Plot 1: Numerical vs Analytical
subplot(4,1,1);
plot(t, euler_sol(1,:),'-ob', 'DisplayName', 'Euler Method');
hold on;
plot(t, heun_sol(1,:),'-+g', 'DisplayName', 'Heun''s Method');
plot(t, rk4_sol(1,:),'-xr', 'DisplayName', 'RK4 Method');
plot(t, analytical_sol, 'k--', 'DisplayName', 'Analytical Solution');
xlabel('Time');
ylabel('Displacement');
title('Numerical Solutions vs Analytical Solution');
legend;
grid on;

% Plot 2: Difference for Euler Method
subplot(4,1,2);
plot(t, euler_diff, 'b');
xlabel('Time');
ylabel('Difference');
title('Difference between Euler Method and Analytical Solution');
grid on;

% Plot 3: Difference for Heun's Method
subplot(4,1,3);
plot(t, heun_diff, 'g');
xlabel('Time');
ylabel('Difference');
title('Difference between Heun''s Method and Analytical Solution');
grid on;

% Plot 4: Difference for RK4 Method
subplot(4,1,4);
plot(t, rk4_diff, 'k');
xlabel('Time');
ylabel('Difference');
title('Difference Between RK4 Method and Analytical Solution');
grid on;

% Define the function that represents the system of ODEs
function dydt = f(t, y, m, k)
    dydt = [y(2); -k/m * y(1)];
end

%%Double Pendulum Analysis

% COMBINED_PENDULUM_ANALYSIS.M
% This script solves the overdriven pendulum problem using Euler, Heun,
% and Runge-Kutta 4th-order (RK4) methods, plots the results, and generates
% a bifurcation diagram using the chosen numerical method.

%% Parameters
alpha = 0.1;   % Damping coefficient
beta = 1.0;    % Gravity term
omega = 2/3;   % Driving frequency
t0 = 0;        % Initial time
tf = 100;      % Final time for bifurcation analysis
dt = 0.01;     % Time step

% Initial conditions
theta_0 = [0.2; 0.0]; % Initial state [theta1; theta2]

% Time vector for the main simulation
t_main = t0:dt:tf;
n_main = length(t_main) - 1;

% Initialize solution vectors
theta_euler = zeros(2, n_main+1);
theta_heun = zeros(2, n_main+1);
theta_rk4 = zeros(2, n_main+1);

% Set initial conditions
theta_euler(:, 1) = theta_0;
theta_heun(:, 1) = theta_0;
theta_rk4(:, 1) = theta_0;

%% Integration using Euler, Heun, and RK4 methods (for initial analysis)
for i = 1:n_main
    theta_euler(:, i+1) = euler_step(theta_euler(:, i), dt, t_main(i), alpha, beta, gamma, omega);
    theta_heun(:, i+1) = heun_step(theta_heun(:, i), dt, t_main(i), alpha, beta, gamma, omega);
    theta_rk4(:, i+1) = rk4_step(theta_rk4(:, i), dt, t_main(i), alpha, beta, gamma, omega);
end

%% Plot results for initial analysis
figure;
hold on;
plot(t_main, theta_euler(1, :), 'r-', 'DisplayName', 'Euler');
plot(t_main, theta_heun(1, :), 'g--', 'DisplayName', 'Heun');
plot(t_main, theta_rk4(1, :), 'b:', 'DisplayName', 'RK4');
legend;
xlabel('Time (s)');
ylabel('Angular Displacement (rad)');
title('Overdriven Pendulum using Euler, Heun, and RK4 Methods');
grid on;
hold off;

%% Bifurcation Diagram using RK4 method (or chosen method)
% Parameters for bifurcation analysis
gamma_start = 1.35;
gamma_end = 1.5;
gamma_step = 0.001;
gamma_values = gamma_start:gamma_step:gamma_end;
num_gammas = length(gamma_values);

% Preallocate for efficiency
theta_samples = cell(num_gammas, 1);

% Time vector for bifurcation analysis
t_bifurcation = t0:dt:tf;
num_steps = length(t_bifurcation);

% Define stroboscopic sampling times (once every driving period after transient)
transient_time = 100; % Time to allow transients to decay
sampling_interval = (2*pi)/omega; % Sampling once every period of the driving force
sampling_times = transient_time:sampling_interval:tf;

% Choose the method for bifurcation analysis
% 1 = Euler, 2 = Heun, 3 = RK4
method_choice = 3;

%% Main Loop: Iterate over gamma values
for idx = 1:num_gammas
    gamma = gamma_values(idx);
    
    % Initialize state vector
    theta = zeros(2, num_steps);
    theta(:,1) = theta_0;
    
    % Integration using the selected method
    for i = 1:num_steps-1
        switch method_choice
            case 1
                theta(:, i+1) = euler_step(theta(:, i), dt, t_bifurcation(i), alpha, beta, gamma, omega);
            case 2
                theta(:, i+1) = heun_step(theta(:, i), dt, t_bifurcation(i), alpha, beta, gamma, omega);
            case 3
                theta(:, i+1) = rk4_step(theta(:, i), dt, t_bifurcation(i), alpha, beta, gamma, omega);
        end
        
        % Keep theta within [-pi, pi] for better visualization
        theta(1, i+1) = mod(theta(1, i+1) + pi, 2*pi) - pi;
    end
    
    % Interpolate theta values at sampling times
    theta_interp = interp1(t_bifurcation, theta(1,:), sampling_times);
    
    % Store samples for current gamma
    theta_samples{idx} = theta_interp;
end

%% Plot Bifurcation Diagram
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
hold on;

for idx = 1:num_gammas
    gamma = gamma_values(idx);
    theta_vals = theta_samples{idx};
    
    % Plot each sample point
    plot(gamma * ones(size(theta_vals)), theta_vals, 'k.', 'MarkerSize', 1);
end

xlabel('\gamma (Driving Force Amplitude)');
ylabel('\theta (Angular Displacement)');
title('Bifurcation Diagram of Overdriven Pendulum');
grid on;
hold off;

%% Function Definitions
function dtheta_dt = pendulum_system(t, theta, alpha, beta, gamma, omega)
    % PENDULUM_SYSTEM Defines the system of ODEs for the overdriven pendulum.
    dtheta_dt = zeros(2,1);
    dtheta_dt(1) = theta(2);
    dtheta_dt(2) = -alpha * theta(2) - beta * sin(theta(1)) + gamma * cos(omega * t);
end

function theta_next = euler_step(theta, dt, t, alpha, beta, gamma, omega)
    % EULER_STEP Performs one step of Euler's method.
    dtheta_dt = pendulum_system(t, theta, alpha, beta, gamma, omega);
    theta_next = theta + dt * dtheta_dt;
end

function theta_next = heun_step(theta, dt, t, alpha, beta, gamma, omega)
    % HEUN_STEP Performs one step of Heun's method (Improved Euler).
    k1 = pendulum_system(t, theta, alpha, beta, gamma, omega);
    k2 = pendulum_system(t + dt, theta + dt * k1, alpha, beta, gamma, omega);
    theta_next = theta + 0.5 * dt * (k1 + k2);
end

function theta_next = rk4_step(theta, dt, t, alpha, beta, gamma, omega)
    % RK4_STEP Performs one step of the 4th-order Runge-Kutta method.
    k1 = pendulum_system(t, theta, alpha, beta, gamma, omega);
    k2 = pendulum_system(t + 0.5 * dt, theta + 0.5 * dt * k1, alpha, beta, gamma, omega);
    k3 = pendulum_system(t + 0.5 * dt, theta + 0.5 * dt * k2, alpha, beta, gamma, omega);
    k4 = pendulum_system(t + dt, theta + dt * k3, alpha, beta, gamma, omega);
    theta_next = theta + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

% Display runtime
toc