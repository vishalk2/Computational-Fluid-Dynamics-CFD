% 1D Steady Heat conduction with Heat Generation & Convection
clear all; clc

% Known Parameters
Lx = 0.04; % Length of x-domian (0.04 m)
M = 21; % Number of grid points chosen along x-domain
T_inf = 30; % Free air temperature (30 C)
h = 45; % Convection Heat transfer coefficient
k = 28; % Thermal conductivity
q = 5e6; % Heat Generation (W/m^3)

dx = Lx/(M-1); % grid size

% Initializing Temperature Matrix
T = zeros(1,M);

% Boundary Conditions
T(1,1) = 0; % Left Boundary Condition

% Computation
err = 1; % Error
epsilon = 1e-8; % Tolerance
n = 0; % Iterator

while err>epsilon
    T_old = T; % T_old = T from kth iteration
    for i=2:M-1
        T(1,i) = (1/2)*(T(1,i+1)+T(1,i-1)+(q/k)*(dx^2));
    end
    % Right Boundary Condition
    T(1,M) = (4*T(1,M-1)-T(1,M-2)+T_inf*(2*h*dx/k))/(3+(2*h*dx/k));
    err = 0; % Error
    for i=1:M
        err = err + power(T(1,i)-T_old(1,i),2);
    end
    err = sqrt(err/M);
    n = n+1; % iterator
end
T
n
err

% Plotting
x = linspace(0,Lx,M);
plot(x,T,'ro-')
xlabel('x-Domain Length'),ylabel('Temperature')
