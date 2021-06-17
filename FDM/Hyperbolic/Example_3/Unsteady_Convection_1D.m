% Example-3 : 1D Linear Convection
% Comparison between FTBS (explicit), Lax-Wendroff (explicit)
% and Lax (explicit) schemes
clear all; clc

% Known Values
Lx = 300; % Length of x-domain
C = 300; % Velocity (300 m/s)
M = 61; % Number of grid points
x = linspace(0,Lx,M);

% Calculated values
dx = Lx/(M-1); % Grid size
Dt = [0.01666,0.015,0.0075]; % Time Steps
CL = zeros(size(Dt)); % Courant Number for each step

% Initial Conditions
% U(x) = 0                      for 0<=x<=50
% U(x) = 100*sin(pi*(x-50)/60)  for 50<=x<=110
% U(x) = 0                      for 110<=x<=300

% Boundary Conditions
% U(x,t) = 0 for x = 0
% U(x,t) = 0 for x = Lx

% Velocity Matrices
U_ftbs = [];
U_lw = [];
U_lax = [];

% Computation
for j=1:length(Dt)
    dt = Dt(j); % Time Step
    lambda = C*dt/dx; % Courant Number
    CL(j) = lambda; % Courant Number for each time step
    
    % Initializing Velocity Matrix
    U_FTBS = zeros(1,M); % U for FTBS
    U_LW = zeros(1,M); % U for Lax-Wendroff
    U_LAX = zeros(1,M); % U for FTCS
    
    % Applying given conditions
    for p=1:M
        if((50<=x(p))&&(x(p)<=110))
            U_FTBS(1,p) = 100*sin(pi*(x(p)-50)/60);
            U_LW(1,p) = 100*sin(pi*(x(p)-50)/60);
            U_LAX(1,p) = 100*sin(pi*(x(p)-50)/60);
        end
    end
    
    % Initial Velocity Profile (remains same for each case)
    U_initial = U_FTBS;
    
    % Computation using each method
    t = 0; % Time
    while t<=0.45
        
        U_FTBS_old = U_FTBS; % FTBS Scheme
        U_LW_old = U_LW; % Lax-Wendroff Scheme
        U_LAX_old = U_LAX; % FTCS Scheme
        
        for i=2:M-1
            U_FTBS(1,i) = (1-lambda)*U_FTBS_old(1,i)+lambda*U_FTBS_old(1,i-1); % FTBS Scheme
            U_LW(1,i) = (1-lambda^2)*U_LW_old(1,i)-(lambda/2)*(1-lambda)*U_LW_old(1,i+1)+(lambda/2)*(1+lambda)*U_LW_old(1,i-1); % Lax-Wendroff Scheme
            U_LAX(1,i) = ((1-lambda)/2)*U_LAX_old(1,i+1)+((1+lambda)/2)*U_LAX_old(1,i-1); % Lax Scheme
        end
        t = t+dt;
    end
    U_ftbs = [U_ftbs;U_FTBS]; % Updating U_ftbs array for each time step chosen
    U_lw = [U_lw;U_LW]; % Updating U_lw array for each time step chosen
    U_lax = [U_lax;U_LAX]; % Updating U_ftcs array for each time step chosen
end

% Plotting

% Results using FTBS
figure;
plot(x,U_initial,'black')
hold on;
plot(x,U_ftbs)
xlabel('Length (x)'),ylabel('Velocity')
title('Velocity Profile')
legend('Exact (C=1)',strcat('C =',num2str(CL(1))),strcat('C =',num2str(CL(2))),strcat('C =',num2str(CL(3))),'Location','bestoutside')
hold off;

% Results using Lax-Wendroff
figure;
plot(x,U_initial,'black')
hold on;
plot(x,U_lw)
xlabel('Length (x)'),ylabel('Velocity')
title('Velocity Profile')
legend('Exact (C=1)',strcat('C =',num2str(CL(1))),strcat('C =',num2str(CL(2))),strcat('C =',num2str(CL(3))),'Location','bestoutside')
hold off;

% Results using FTCS
figure;
plot(x,U_initial,'black')
hold on;
plot(x,U_lax)
xlabel('Length (x)'),ylabel('Velocity')
title('Velocity Profile')
legend('Exact (C=1)',strcat('C =',num2str(CL(1))),strcat('C =',num2str(CL(2))),strcat('C =',num2str(CL(3))),'Location','bestoutside')
hold off;
