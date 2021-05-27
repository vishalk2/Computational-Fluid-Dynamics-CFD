% Steady 1D Heat Conduction Example-1
clear all; clc

% Known Values
Lx = 1; % Length of x-domain
M = 21; % Number of Points

% Calculated values
dx = Lx/(M-1); % Step Size

% Initilizing Temperature Matrix
T = zeros(1,M);

% Boundary Conditions
T(1,1) = 100; % Left Boundary T = 100 C
T(1,M) = 0; % Right Boundary T = 0 C

% Computation
epsilon = 1e-8; % Error Sensivity
err = 1; % Error
n = 0; % Iterator

while err>epsilon
    T_old = T;
    for i=2:M-1
        T(1,i) = (1/2)*(T(1,i+1)+T(1,i-1));
    end
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
plot(x,T,'r-o')
xlabel('Domain Length (x)'),ylabel('Temperature')
title('Temperature Variation')
