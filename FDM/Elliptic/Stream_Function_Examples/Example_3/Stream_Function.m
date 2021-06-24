% Example-3 : Stream Function Variation
% Jacobi VS Gauss-Seidel VS PSOR
clear all; clc

% Known Values
Lx = 40; % Length of x-domain
Ly = 40; % Length of y-domain
M = 81; % Number of grid points in x-domain
N = 81; % Number of grid points in y-domain
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);

% Calculated values
dx = Lx/(M-1);
dy = Ly/(N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Jacobi Method
% Coefficients
ap = 2*((1/(dx^2))+(1/(dy^2)));
ae = 1/(dx^2);
aw = 1/(dx^2);
an = 1/(dy^2);
as = 1/(dy^2);

% Initializing Temperature Matrix
Psi_Jacobi = zeros(N,M);
Psi_old_Jacobi = zeros(N,M);

% Boundary Conditions
Psi_Jacobi(N,:) = 0; % Bottom Boundary
Psi_Jacobi(1,:) = 0; % Top Boundary
Psi_Jacobi(:,M) = 0; % Right Boundary

% Computation - Jacobi Method
epsilon = 1e-8; % Tolerance or Error Sensitivity
err_Jacobi = 1; % Error
Err_Jacobi = []; % Error Array
n_Jacobi = 0; % Iterator

while err_Jacobi>epsilon
    % Left Boundary Condition
    for j=1:N
        if(10<=y(j) && y(j)<=30)
            Psi_Jacobi(j,1) = 1;
        else
            Psi_Jacobi(j,1) = 0;
        end
    end
    
    for j=1:N
        for i=1:M
            Psi_old_Jacobi(j,i) = Psi_Jacobi(j,i);
        end
    end
    for j=2:N-1
        for i=2:M-1
            Psi_Jacobi(j,i) = (1/ap)*(ae*Psi_old_Jacobi(j+1,i) + aw*Psi_old_Jacobi(j-1,i) + an*Psi_old_Jacobi(j,i+1) + as*Psi_old_Jacobi(j,i-1));
        end
    end
    err_Jacobi = 0; % Error
    for j=1:N
        for i=1:M
            err_Jacobi = err_Jacobi + power(Psi_Jacobi(j,i)-Psi_old_Jacobi(j,i),2);
        end
    end
    err_Jacobi = sqrt(err_Jacobi/(M*N));
    Err_Jacobi(end+1) = err_Jacobi;
    n_Jacobi = n_Jacobi+1; % Iterator
end

% Plotting using Jacobi Method Results
figure;
[X,Y] = meshgrid(x,y);
contourf(X,Y,Psi_Jacobi,':','ShowText','on'),colorbar,colormap(jet)
xlabel('X'),ylabel('Y')
title({['Stream Function Contour Plot'];['Jacobi Method'];['Iterations taken : ',num2str(n_Jacobi)]})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. Gauss-Seidel Method
% Initializing Temperature Matrix
Psi_GS = zeros(N,M);

% Boundary Conditions
Psi_GS(N,:) = 0; % Bottom Boundary
Psi_GS(1,:) = 0; % Top Boundary
Psi_GS(:,M) = 0; % Right Boundary

% Computation - Gauss Seidel Method
epsilon = 1e-8; % Tolerance or Error Sensitivity
err_GS = 1; % Error
Err_GS = []; % Error Array
n_GS = 0; % Iterator

while err_GS>epsilon
    % Left Boundary Condition
    for j=1:N
        if(10<=y(j) && y(j)<=30)
            Psi_GS(j,1) = 1;
        else
            Psi_GS(j,1) = 0;
        end
    end
    
    err_GS = 0;
    for j=2:N-1
        for i=2:M-1
            Psi_GS_temp = Psi_GS(j,i);
            Psi_GS(j,i) = (1/ap)*(ae*Psi_GS(j+1,i) + aw*Psi_GS(j-1,i) + an*Psi_GS(j,i+1) + as*Psi_GS(j,i-1));
            err_GS = err_GS + power(Psi_GS(j,i)-Psi_GS_temp,2);
        end
    end
    err_GS = sqrt(err_GS/(M*N));
    Err_GS(end+1) = err_GS;
    n_GS = n_GS+1; % Iterator
end

% Plotting using Gauss-Seidel Method Results
figure;
[X,Y] = meshgrid(x,y);
contourf(X,Y,Psi_GS,':','ShowText','on'),colorbar,colormap(jet)
xlabel('X'),ylabel('Y')
title({['Stream Function Contour Plot'];['Gauss-Seidel Method'];['Iterations taken : ',num2str(n_GS)]})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. PSOR Method
% Initializing Temperature Matrix
Psi_PSOR = zeros(N,M);
Psi_old_PSOR = zeros(N,M);

% Boundary Conditions
Psi_PSOR(N,:) = 0; % Bottom Boundary
Psi_PSOR(1,:) = 0; % Top Boundary
Psi_PSOR(:,M) = 0; % Right Boundary

% Finding Optimum Relaxation Factor
B = dx/dy;
a = ((cos(pi/(M-1))+(B^2)*cos(pi/(N-1)))/(1+(B^2)))^2;
w = (2-sqrt(1-a))/a; % Optimum Relaxation Factor

% Computation
epsilon = 1e-8; % Error Sensivity
err_PSOR = 1; % Error
Err_PSOR = []; % Error Array
n_PSOR = 0; % Iterator

while err_PSOR>epsilon
    % Left Boundary Condition
    for j=1:N
        if(10<=y(j) && y(j)<=30)
            Psi_PSOR(j,1) = 1;
        else
            Psi_PSOR(j,1) = 0;
        end
    end
    
    for j=1:N
        for i=1:M
            Psi_old_PSOR(j,i) = Psi_PSOR(j,i);
        end
    end
    for j=2:N-1 % For Internal Grid points
        for i=2:M-1 % For Internal Grid points
            % T(j,i)_(k+1)th_iteration Calculation
            Psi_PSOR(j,i) = (1-w)*Psi_old_PSOR(j,i)+(w/(2*(1+B^2)))*((B^2)*Psi_PSOR(j,i-1)+Psi_PSOR(j-1,i)+Psi_old_PSOR(j,i+1)+(B^2)*Psi_old_PSOR(j+1,i));
        end
    end
    err_PSOR = 0; % Error
    for j=1:N
        for i=1:M
            err_PSOR = err_PSOR + power(Psi_PSOR(j,i)-Psi_old_PSOR(j,i),2);
        end
    end
    err_PSOR = sqrt(err_PSOR/(M*N));
    Err_PSOR(end+1) = err_PSOR;
    n_PSOR = n_PSOR+1; % Iterator
end

% Plotting using PSOR Method Results
figure;
[X,Y] = meshgrid(x,y);
contourf(X,Y,Psi_PSOR,':','ShowText','on'),colorbar,colormap(jet)
xlabel('X'),ylabel('Y')
title({['Stream Function Contour Plot'];['PSOR Method'];['Iterations taken : ',num2str(n_PSOR)]})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Errors VS Iterations Comparison Plot
figure;
x_Jacobi = 1:1:n_Jacobi;
x_GS = 1:1:n_GS;
x_PSOR = 1:1:n_PSOR;
loglog(x_Jacobi,Err_Jacobi,'r-',x_GS,Err_GS,'b-',x_PSOR,Err_PSOR,'g-')
hold on;
n = 1:1:max([n_Jacobi,n_GS,n_PSOR]);
Tolerance = (1e-8)*ones(1,max([n_Jacobi,n_GS,n_PSOR]));
plot(n,Tolerance,'k-','LineWidth',2);
ylim([1e-10,1e0])
hold off;
xlabel('log(Iterations)'),ylabel('log(Error)')
title('Error VS Iterations Comparison between Schemes')
legend('Jacobi','Gauss-Seidel','PSOR','Tolerance')
