% Example-2 : 1D Linear Convection
clear all; clc

% Velocity is 2m/s in the range x=0.1 to x=0.3 and is
% 1m/s everywhere else.
Xs = 0.1; % Starting point of square wave
Xe = 0.3; % Ending point of square wave

% Simulating based on different grid sizes

% Known Values
Lx = 1; % Length of x-domain
C = 1; % Velocity (1 m/s)
M = [21,41,81,161]; % Number of Grid Points
Err = zeros(size(M)); % Error from each grid size

for j=1:length(M)
    % Calculated Values
    dx = Lx/(M(j)-1); % Grid Size for each "M"
    dt = 0.01; % Time Step
    x = linspace(0,Lx,M(j)); % X-Domain
    lambda = C*dt/dx; % Courant Number
    
    % Initializing Velocity matrix
    U = ones(1,M(j));
    
    % Applying Given conditions
    Ms(j) = round((Xs/dx)+1); % Starting node for square wave
    Me(j) = round((Xe/dx)+1); % Ending node for square qave
    U(:,Ms(j):Me(j)) = 2; % Velocity of wave in x=0.1 to x=0.3 range
    
    U_initial = U; % Initial Velocity
    
    % Computation
    t = 0; % Time
    n = 0; % Iterator
    while t<=0.4
        err = 0; % Error
        U_old = U;
        for i=2:M(j)
            U(1,i) = (1-lambda)*U_old(1,i)+lambda*U_old(1,i-1);
            err = err + power(U(1,i)-U_old(1,i),2); % Error
        end
        err = sqrt(err/M(j)); % Error calculation
        t = t+dt;
        n = n+1; % Iterator
        
        % Plotting
        figure(j);
        plot(x,U_initial,'r-')
        %axis([0 1 0 3])
        hold on
        plot(x,U,'b-')
        xlabel('Length (x)'),ylabel('Velocity (U)')
        title('Velocity Profile')
        legend('Exact Numerical Velocity Profile','Numerical Velocity Profile','Location','northeastoutside')
    end
    Err(j) = err; % Error for each gid size variation
end
