% Example-4 : 2nd Order Wave Equation - 1D
% Solving using Mid-Point Leap Frog Method
clear all; clc

% Known Values
Lx = 300; % Length of x-domain
a = 300; % Velocity (in m/s)
M = 61; % Number of grid points
x = linspace(0,Lx,M);

% Calculated values
dx = Lx/(M-1); % Grid size
Dt = [0.01666,0.015,0.0075]; % Time Steps
CL = zeros(size(Dt)); % Courant Number for each step

% Initial Conditions (@ t=0)
% U(x) = 0                        for 0<=x<=20
% U(x) = 100*sin(pi*(x-20)/120)   for 20<=x<=140
% U(x) = 0                        for 140<=x<=180
% U(x) = 100*sin(pi*(x-180)/120)  for 180<=x<=300
% dU/dt = 0

% Boundary Conditions
% U(x,t) = 0 for x = 0
% U(x,t) = 0 for x = Lx

% Velocity Matrices
U_array = [];

% Computation
for j=1:length(Dt)
    dt = Dt(j); % Time Step
    lambda = a*dt/dx; % Courant Number
    CL(j) = lambda; % Courant Number for each time step
    
    % Initializing Velocity Matrix
    U = zeros(1,M); % U for 
    
    % Applying given conditions
    for p=1:M
        if((20<=x(p))&&(x(p)<=140))
            U(1,p) = 100*sin(pi*(x(p)-20)/120);
        end
        if((180<=x(p))&&(x(p)<=300))
            U(1,p) = 100*sin(pi*(x(p)-180)/120);
        end
    end
    
    % Initial Velocity Profile
    U_initial = U;
    
    % Computation using Mid-Point Leap Frog Method
    t = 0; % Time
    while t<=0.28
        U_old = U;
        for i=2:M-1
            U(1,i) = U_old(1,i)+((lambda^2)/2)*(U_old(1,i-1)-2*U_old(1,i)+U_old(1,i+1));
        end
        t = t+dt;
    end
    U_array = [U_array;U]; % Updating U_array for each time step chosen
end

% Plotting
figure;
plot(x,U_initial,'k')
hold on;
grid on;
plot(x,U_array)
ylim([-10 110])
xlabel('Length (x)'),ylabel('Velocity (U)')
title('Velocity Profile','Time = 0.28')
set(gca,'xtick',[0:50:Lx])
set(gca,'ytick',[0:20:100])
legend('Exact (C=1)',strcat('C =',num2str(CL(1))),strcat('C =',num2str(CL(2))),strcat('C =',num2str(CL(3))),'Location','bestoutside')
hold off;
