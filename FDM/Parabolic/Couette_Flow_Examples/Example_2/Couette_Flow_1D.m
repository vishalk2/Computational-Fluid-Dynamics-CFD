% 1D Couette Flow Example-2 using Crank-Nicolson Scheme
% Fixed Plate -> Top Plate
% Moving Plate -> Bottom Plate
clear all; clc

% Known parameters
h = 0.04; % Distance between the plates (0.04 m)
N = 41; % No. of grid points considered
dy = h/(N-1); % Grid Size
U_m = 40; % Velocity of Moving Plate (40 m/s)
v = 0.000217; % v = Uh/Re; Re = reynolds number (m^2/s)
dt = 2e-3; % Time Step
gamma_y = v*dt/(dy^2); % Diffusion Factor Coefficient

% Initialization
t = 0; % Time
u = zeros(N,1); % Velocity Matrix
U = []; % Velocities at particular time steps
T = []; % Particular Time Steps
n = 0; % Iterator

% Initial & Boundary Conditions
u(N,1) = U_m; % Moving Plate (Bottom Plate)
u(1,1) = 0; % Stationary Plate (Top Plate)

% Computation - Implicit Method
% Crank Nicolson Method
while t<=10
    u_old = u; % u_old at nth time step & u at (n+1)th time step
    for i=2:N-1
        u(i,1) = ((gamma_y/(2+2*gamma_y))*(u(i+1,1)+u(i-1,1)+u_old(i+1,1)+u_old(i-1,1)))+((1-gamma_y)/(1+gamma_y))*u_old(i,1);
    end
    U = [U;u'];
    T(end+1) = t;
    t = t+dt; % Time
    n = n+1; % Iterator
end
u
n
t

% Plotting
figure;
hold on;
grid on;
%plot(u,h:-dy:0)
for i=1:150:n
    plot(U(i,:),h:-dy:0,'b-')
end
xlabel('Velocity (u)'),ylabel('Height (Y)')
title('Flow between 2 parallel plates')
hold off;
