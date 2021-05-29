% Example-1
% Couette Flow : Flow between 2 infinitely parallel plates
clear all; clc

% Known Parameters
N = 101; % No. of grid points
h = 1; % Distance between the plates (Length of y-domain)
dy = h/(N-1); % Grid Size
u_m = 1; % Velocity of Moving plate
Re = 100; % Reynold's Number (Aribitrary)

% Calculations
% Since nu = gamma = U*h/Re
gamma = (u_m*h/Re); % Diffusion Coefficient

% Condition for stability: gamma_y <= (1/2)
% gamma_y = nu*dt/(dy^2)
dt = (dy^2)/(2*gamma); % Time Step
gamma_y = gamma*dt/(dy*dy);

% Initialization
u = zeros(N,1); % Velocity Matrix
U = []; % Velocity variation array with time step
t = 0; % Time
T = []; % Time array

% Initial & Boundary Conditions
u(1,1) = u_m; % Moving plate
u(N,1) = 0; % Stationary plate

% Computation - Explicit Method
% Forward in Time Central in Space (FTCS)
n = 0; % Iterator
while t<=50
    u_old = u;
    for j=2:N-1 % For Internal Nodes
        u(j,1) = u_old(j,1) + gamma_y*(u_old(j+1,1)-2*u_old(j,1)+u_old(j-1,1));
    end
    n = n+1;
    t = t + dt;
    U = [U;u'];
    T(end+1) = t;
end
u
n
t

% Plotting
figure;
hold on;
grid on;
for i=1:(2/dt):n
    plot(U(i,:),h:-dy:0,'b-')
end
xlabel('Velocity (u)'),ylabel('Height (Y)')
title('Flow between 2 parallel plates')
hold off;
