% 1D heat Conduction Example-1 : Fin connected to a wall
% Tip of fin is insulated
clear all; clc

% Known Parameters
Lx = 1; % Length of x-domain (1 metre)
k = 1.2; % Thermal Conductivity of wall
h = 8; % Convection heat Transfer Coefficient
d = 0.7; % Diameter of Fin
M = 21; % Number of Points

% Calculations
dx = Lx/(M-1);
P = (22/7)*d; % Perimeter
A = (22/28)*(d^2); % Cross sectional area of Fin
m = sqrt(h*P/(k*A));
T_inf = 300; % Free Air Temperature; T_infinity = 300K

% Initializing temperature Matrix
T = zeros(1,M);

% Boundary Conditions
T(1,1) = 500; % Left Boundary; T = 500K

err = 1;
epsilon = 1e-8;
n = 0; % Iterator

while err>epsilon
    T_old = T;
    for i=2:M-1
        T(1,i) = (1/(2+(m^2)*(dx^2)))*(T(1,i+1)+T(1,i-1)+(m^2)*(dx^2)*T_inf);
    end
    T(1,M) = (4*T(1,M-1)-T(1,M-2))/3; % Right Adiabatic Boundary Condition using 2nd order approximation
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
plot(x,T,'r-o');
xlabel('x-Domain Length'),ylabel('Temperature')
title('Temperature Variation')
