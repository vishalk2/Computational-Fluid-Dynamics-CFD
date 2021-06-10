% Example-7 : 2D Heat Conduction & Convection using PSOR
clear all; clc

% Known Values
Lx = 0.2; % Length of x-domain
Ly = 0.1; % Length of y-domain
M = 41; % Number of Grid Points along x-direction
N = 21; % Number of Grid Points along y-direction
hg = 1000; % Heat Convection Coefficient of gases
Tg = 2000; % Temperature of gases
hb = 8000; % Heat Convection Coefficient of coolant
Tb = 60; % Temperature of coolant
k = 20; % Thermal Conductivity

% Calculated Values
dx = Lx/(M-1); % Grid size in x-direction
dy = Ly/(N-1); % Grid size in y-direction
B = dx/dy;
%w = 1.9; % Relaxation Factor
W = 1:0.001:2;
NI = zeros(size(W));

for p=1:length(W)
    w = W(p);
    
    % Computation
    epsilon = 1e-8; % Error Sensivity
    err = 1; % Error
    n = 0; % Iterator
    
    % Initializing Temperature Matrix
    T = zeros(N,M);
    
    while err>epsilon
        % Applying Boundary Conditions
        for i=1:M % Top Boundary Condition (Convective)
            T(1,i) = (4*T(2,i)-T(3,i)+(2*dy*hg/k)*Tg)/(3+(2*dy*hg/k));
        end
        for j=1:N
            T(j,1) = (4*T(j,2)-T(j,3))/3; % Left Boundary Condition (Adiabatic)
            T(j,M) = (4*T(j,M-1)-T(j,M-2))/3; % Right Boundary Condition (Adiabatic)
        end
        for i=1:21
            T(N,i) = (4*T(N-1,i)-T(N-2,i))/3; % Half Bottom Boundary Condition (Adiabatic)
        end
        for i=21:M
            T(N,i) = (4*T(N-1,i)-T(N-2,i)+(2*dy*hb/k)*Tb)/(3+(2*dy*hb/k));
        end
        
        % Updating T_old i.e. T from previous iteration
        for j=1:N
            for i=1:M
                T_old(j,i) = T(j,i);
            end
        end
        % Updating T from current iteration
        for j=2:N-1
            for i=2:M-1
                T(j,i) = (1-w)*T_old(j,i)+(w/(2*(1+B^2)))*((B^2)*T(j,i-1)+T(j-1,i)+T_old(j,i+1)+(B^2)*T_old(j+1,i));
            end
        end
        err = 0; % Error
        for j=1:N
            for i=1:M
                err = err + power(T(j,i)-T_old(j,i),2);
            end
        end
        err = sqrt(err/(M*N));
        n = n+1; % Iterator
    end
    NI(p) = n;
end

% Plotting
plot(W,NI,'b-'),xlim([0.8,2.2])
xlabel('Relaxation Factor'),ylabel('Number of Iterations')
title('Number of Iterations taken VS Relaxation factor')
hold on;
plot(W(islocalmin(NI)),NI(islocalmin(NI)),'ro')
hold off;
