% Example-5 : Non-Linear Wave Equation - 1D
% Lax Method
clear all; clc

% Known Values
Lx = 4; % Length of x-domain
M = 81; % Number of grid points
x = linspace(0,Lx,M); % x-Domain
U_max = 1; % Max Velocity

% Calculated Values
dx = Lx/(M-1); % Grid Size in x-direction
Dt = [dx,dx/2]; % Time step chosen so in order to obey |(dt/dx)*U_max| <= 1 for stability
T = [0.5,1.0,1.5,2.0]; % Time Intervals chosen

% Computation
for j=1:length(Dt)
    dt = Dt(j); % Time Step
    N = T./dt; % Iteration Number at each time interval
    
    % Initializing Velocity Matrix
    U = zeros(1,M);
    % Velocity Matrices
    U_array = [];
    
    % Applying given initial & boundary conditions
    for p=1:M
        if((0<=x(p))&&(x(p)<=2))
            U(1,p) = 1;
        end
        if((2<=x(p))&&(x(p)<=4))
            U(1,p) = 0;
        end
    end
    
    % Initial Velocity Profile
    U_initial = U;
    
    % Computation using Lax Method (FTCS)
    t = 0; % Time
    n = 0; % Iterator
    while t<=2
        U_old = U;
        for i=2:M-1
            U(1,i) = ((U_old(1,i+1)+U_old(1,i-1))/2)-(dt/(4*dx))*((U_old(1,i+1)^2)-(U_old(1,i-1)^2));
        end
        t = t+dt;
        n = n+1; % Iterator
        U_array = [U_array;U]; % Updating U_array
    end
    
    % Plotting
    figure;
    plot(x,U_initial,'k',"LineWidth",2)
    hold on;
    grid on;
    ylim([-0.5 1.5])
    xlabel('Length (x)'),ylabel('Velocity (U)')
    title('Velocity Profile - Lax Method',strcat('Time Step :',num2str(dt)))
    set(gca,'xtick',[0:1:Lx])
    set(gca,'ytick',[-0.5:0.5:1.5]) % U_max = 1
    
    for k=1:n
        for p=1:length(N)
            if k==N(p)
                plot(x,U_array(k,:),'color',[0.9100 0.4100 0.1700])
            end
        end
    end
    legend('T = 0.0 sec',strcat('T =',num2str(T(1))),strcat('T =',num2str(T(2))),strcat('T =',num2str(T(3))),strcat('T =',num2str(T(4))),'Location','bestoutside')
    hold off;
end
