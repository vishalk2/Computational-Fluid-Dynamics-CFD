% 1D Couette Flow Example-2
% Fixed Plate -> Top Plate
% Moving Plate -> Bottom Plate
clear all; clc

% Known parameters
h = 0.04; % Distance between the plates
N = 41; % No. of grid points considered
dy = h/(N-1); % Grid Size
U_m = 40; % Velocity of Moving Plate
v = 0.000217; % v = 1/Re; Re = reynolds number
dt = 0.002; % Time Step
gamma_y = v*dt/(dy^2); % Diffusion Factor Coefficient

% Initialization
t = 0; % Time
u = zeros(N,1); % Velocity Matrix
U = []; % Velocities at particular time steps
T = []; % Particular Time Steps

% Initial & Boundary Conditions
u(N,1) = U_m; % Moving Plate
u(1,1) = 0; % Stationary Plate

% Computation - Implicit Method
% Backward in Time, Central in Space (BTCS)
while t<=1
    u_old = u; % u_old at nth time step & u at (n+1)th time step
    for i=2:N-1
        u(i,1) = (1/(2*gamma_y))*(gamma_y*u(i+1,1)+u(i-1,1)+u_old(i,1));
    end
    t = t+dt; % Time
end
