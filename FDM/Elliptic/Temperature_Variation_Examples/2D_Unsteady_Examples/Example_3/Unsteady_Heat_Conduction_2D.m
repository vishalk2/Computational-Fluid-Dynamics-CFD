% Unsteady Heat Conduction 2D
clear all; clc

% Known Values
Lx = 1; % Length of x-domain
Ly = 1; % Length of y-domain
alpha = 1.4; % Thermal Diffusivity
M = 41; % Number of grid points in x-direction
N = 41; % Number of grid points in y-direction

% Calculated Values
dx = Lx/(M-1); % Grid size in x-direction
dy = Ly/(N-1); % Grid size in y-direction
dt = 1e-3; % Time Step
gamma_x = alpha*dx/(dt^2);
gamma_y = alpha*dy/(dt^2);

% Initializing Temperature matrix
T = zeros(N,M);

% Boundary Conditions
T(:,1) = 500; % Left Wall Boundary Condition
T(:,M) = 1000; % Right Wall Boundary Condition

% Computation
n = 0; % Iterator
t = 0;
while t<50
    for i=1:M
        T(1,i) = T(2,i); % Top Neumann Boundary Condition
        T(N,i) = T(N-1,i); % Bottom Neumann Boundary Condition
    end
    err = 0; % Error
    for j=2:N-1
        for i=2:M-1
            T_temp = T(j,i);
            T(j,i) = (1/(1+2*gamma_x+2*gamma_y))*(gamma_x*(T(j+1,i)+T(j-1,i))+gamma_y*(T(j,i+1)+T(j,i-1))+T_temp);
            err = err + power(T(j,i)-T_temp,2);
        end
    end
    err = sqrt(err/(M*N));
    n = n+1;
    t = t+dt;
end
T
n
err

% Plotting
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,T,':','ShowText','on'),colorbar,colormap(jet)
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y'),title('Temperature Variation (T)')
