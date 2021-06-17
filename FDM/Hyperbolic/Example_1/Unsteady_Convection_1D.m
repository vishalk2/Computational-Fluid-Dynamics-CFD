% Example-1 : 1D Linear Convection
clear all; clc

% Velocity is 2m/s in the range x=0.1 to x=0.3 and is
% 1m/s everywhere else.

% Simulating based on different size of time steps

% Known Values
Lx = 1; % Length of x-domain
C = 1; % Velocity (1 m/s)
M = 41; % Number of Grid Points
x = linspace(0,Lx,M);

% Calculated Values
dx = Lx/(M-1); % Grid Size in x-direction
Dt = [1e-1,1e-2,1e-3,1e-4]; % Time Steps
NI = zeros(size(Dt)); % Number of iterations taken for each time step
Err = zeros(size(Dt)); % Error from each time step;

for j=1:length(Dt)
    dt = Dt(j); % Time Step
    lambda = C*dt/dx; % Courant Number
    
    % Initializing Velocity Matrix
    U = ones(1,M);
    
    % Applying given conditions
    U(1,5:13) = 2;
    U(1,1) = 1;
    
    % Computation - Explicit Method (First Upwind Difference Method)
    U_initial = U;
    
    t = 0; % Time
    n = 0; % Iterator
    while t<=0.4
        err = 0; % Error
        U_old = U;
        for i=2:M
            U(1,i) = (1-lambda)*U_old(1,i)+lambda*U_old(1,i-1);
            err = err + power(U(1,i)-U_old(1,i),2); % Error
        end
        err = sqrt(err/M); % Error calculation
        t = t+dt;
        n = n+1; % Iterator
    end
    NI(j) = n;
    Err(j) = err;
    
    % Plotting
    plot(x,U,'LineWidth',1)
    hold on;
    axis([0 1 1 2])
end

% Plotting
plot(x,U_initial,'black','LineWidth',2)
legend('Time Step : 0.1','Time Step : 0.01','Time Step : 0.001','Time Step : 0.0001','Original Profile','Location',"northeastoutside")
hold off;
