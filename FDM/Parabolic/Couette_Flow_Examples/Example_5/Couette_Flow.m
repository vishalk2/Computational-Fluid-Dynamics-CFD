% Example-5 : 1D Flow Variation  using Crank-Nicolson Scheme
clear all; clc

% Known parameters
h = 1; % Distance between the plates
N = 101; % No. of grid points considered
dy = h/(N-1); % Grid Size
U_m = 0; % Velocity of Moving Plate
v = 0.01; % v = Uh/Re; Re = reynolds number (m^2/s)
dt = 5e-3; % Time Step
gamma_y = v*dt/(dy^2); % Diffusion Factor Coefficient
y = 0:dy:h; % Y-domain

% Initialization
t = 0; % Time
u = ones(N,1); % Velocity Matrix
U = []; % Velocities at particular time steps
T = []; % Particular Time Steps
n = 0; % Iterator

% Initial Condition
for j=1:N
    if(0<=y(j) && y(j)<=(1/2))
        u(j,1) = 2*y(j);
    elseif((1/2)<=y(j) && y(j)<=1)
        u(j,1) = 2*(1-y(j));
    end
end

% Computation - Implicit Method
% Crank Nicolson Method
while t<=20
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
%figure;
%hold on;
%grid on;
%plot(u,h:-dy:0)
%for i=1:50:n
%    plot(U(i,:),h:-dy:0,'b-')
%end
%xlabel('Velocity (u)'),ylabel('Height (Y)')
%title('Flow between 2 parallel plates')
%hold off;

% Animation
figure;
myVideo = VideoWriter('CouetteFlow'); %open video file
myVideo.FrameRate = 10;
open(myVideo)

y = h:-dy:0;
for i=1:25:n
    plot(U(i,:),y,'b-')
%    hold on;
    axis([0 1 0 h])
    pause(0.2)
    xlabel('U (in m/s)'),ylabel('Y (in m)')
    frame = getframe(gcf);
    writeVideo(myVideo,frame);
end
close(myVideo)
