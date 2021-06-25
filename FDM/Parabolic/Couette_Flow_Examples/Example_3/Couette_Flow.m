% 1D Couette Flow Example-3 using Crank-Nicolson Scheme
% Fixed Plate -> Bottom Plate
% Moving Plate -> Top Plate
clear all; clc

% Known parameters
h = 5; % Distance between the plates (5 m)
N = 51; % No. of grid points considered
dy = h/(N-1); % Grid Size
U_m = 100; % Velocity of Moving Plate (100 m/s)
v = 1; % v = Uh/Re; Re = reynolds number (m^2/s)
dt = 5e-3; % Time Step
gamma_y = v*dt/(dy^2); % Diffusion Factor Coefficient

% Initialization
t = 0; % Time
u = ones(N,1); % Velocity Matrix
U = []; % Velocities at particular time steps
T = []; % Particular Time Steps
n = 0; % Iterator

% Initial & Boundary Conditions
u = 20*u; % Initial Condition

% Computation - Implicit Method
% Crank Nicolson Method
while t<=10
    if t>0
        u(N,1) = 0; % Stationary Plate Boundary Condition
        u(1,1) = U_m; % Moving Plate Boundary Condition
        u_old = u; % u_old at nth time step & u at (n+1)th time step
        for i=2:N-1
            u(i,1) = ((gamma_y/(2+2*gamma_y))*(u(i+1,1)+u(i-1,1)+u_old(i+1,1)+u_old(i-1,1)))+((1-gamma_y)/(1+gamma_y))*u_old(i,1);
        end
    end
    U = [U;u'];
    T(end+1) = t;
    t = t+dt; % Time
    n = n+1; % Iterator
end
%u
n
t

% Plotting
figure;
hold on;
grid on;
plot(u,h:-dy:0)
for i=1:50:n
    plot(U(i,:),h:-dy:0,'b-')
end
xlabel('Velocity (u)'),ylabel('Height (Y)')
title('Flow between 2 parallel plates')
hold off;

% Animation
figure;
myVideo = VideoWriter('CouetteFlow'); %open video file
myVideo.FrameRate = 10;
open(myVideo)

y = h:-dy:0;
for i=1:25:n
    plot(U(i,:),y,'b-')
%    hold on;
    axis([0 U_m 0 h])
    pause(0.2)
    xlabel('U (in m/s)'),ylabel('Y (in m)')
    frame = getframe(gcf);
    writeVideo(myVideo,frame);
end
close(myVideo)
