% 2D Heat Conduction & Convection Example-1
clear all; clc

% Known Parameters
Lx = 1; % Length of X-Domain
Ly = 1; % Length of Y-Domain
M = 21; % Number of grid points starting from 1 in x-direction
N = 21; % Number of grid points starting from 1 in y-direction
h = 10; % Convective Heat Transfer Coefficient (W/K/m^2)
T_inf = 300; % Free air Temperature T_infi=300K
k = 0.8; % Thermal Conductivity of the wall (W/K/m)

% Calculated constants
dx = Lx/(M-1); % Grid size along x-direction
dy = Ly/(N-1); % Grid size along y-direction
B = dx/dy;

% Coefficients
ap = 2*((1/(dx^2))+(1/(dy^2)));
ae = 1/(dx^2);
aw = 1/(dx^2);
an = 1/(dy^2);
as = 1/(dy^2);

% Initializing Temperature matrix
T = zeros(N,M);

% Boundary Conditions
T(:,1) = 500; % Left Boundary T=500K
T(N,:) = 500; % Bottom Boundary T=500K
T(:,M) = 500; % Right Boundary T=500K

% Computation
epsilon = 1e-8; % Error Sensivity
err = 1; % Error
n = 0; % Iterator

while err>epsilon
    err = 0; % Error
    for j=2:N-1 % For Internal Grid points
        for i=2:M-1 % For Internal Grid points
            T_old = T(j,i);
            T(j,i) = (1/ap)*(ae*T(j+1,i) + aw*T(j-1,i) + an*T(j,i+1) + as*T(j,i-1));
            err = err + power(T(j,i)-T_old,2);
        end
    end
    for i=2:M-1 % Top Neumann Boundary Condition
        T(1,i) = (4*T(2,i)-T(3,i)+(2*dy*(h/k))*T_inf)/(3+(2*dy*(h/k)));
    end
    err = sqrt(err/(M*N));
    n = n+1; % Iterator
end
T
n
err

% Plotting
x = linspace(0,Lx,M);
y = linspace(0,Ly,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,T,'ShowText','on'),colorbar
xlabel('X'),ylabel('Y'),title('Temperature Variation (T)')
