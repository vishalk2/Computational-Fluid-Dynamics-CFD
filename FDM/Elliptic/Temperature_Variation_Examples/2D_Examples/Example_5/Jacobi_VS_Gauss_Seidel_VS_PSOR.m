% Temperature Variation Example-5
% Comparing the Iterations & Results of Jacobi, Gauss-Seidel & PSOR methods

% Known Parameters
Lx = 1; % Length of domain along x-direction
Ly = 1; % Length of domain along y-direction
M = 41; % Number of grid points starting from 1 in x-direction
N = 41; % Number of grid points starting from 1 in y-direction

% Calculated constants
dx = Lx/(M-1); % Grid size along x-direction
dy = Ly/(N-1); % Grid size along y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Jacobi Method
% Coefficients
ap = 2*((1/(dx^2))+(1/(dy^2)));
ae = 1/(dx^2);
aw = 1/(dx^2);
an = 1/(dy^2);
as = 1/(dy^2);

% Initializing Temperature Matrix
T_Jacobi = zeros(N,M);
T_old_Jacobi = zeros(N,M);

% Boundary Conditions
T_Jacobi(N,:) = 900; % Bottom Boundary
T_Jacobi(1,:) = 600; % Top Boundary
T_Jacobi(:,1) = 400; % Left Boundary
T_Jacobi(:,M) = 800; % Right Boundary

% Computation - Jacobi Method
epsilon = 1e-8; % Tolerance or Error Sensitivity
err_Jacobi = 1; % Error
Err_Jacobi = []; % Error Array
n_Jacobi = 0; % Iterator

while err_Jacobi>epsilon
    for j=1:N
        for i=1:M
            T_old_Jacobi(j,i) = T_Jacobi(j,i);
        end
    end
    for j=2:N-1
        for i=2:M-1
            T_Jacobi(j,i) = (1/ap)*(ae*T_old_Jacobi(j+1,i) + aw*T_old_Jacobi(j-1,i) + an*T_old_Jacobi(j,i+1) + as*T_old_Jacobi(j,i-1));
        end
    end
    err_Jacobi = 0; % Error
    for j=1:N
        for i=1:M
            err_Jacobi = err_Jacobi + power(T_Jacobi(j,i)-T_old_Jacobi(j,i),2);
        end
    end
    err_Jacobi = sqrt(err_Jacobi/(M*N));
    Err_Jacobi(end+1) = err_Jacobi;
    n_Jacobi = n_Jacobi+1; % Iterator    
end

% Plotting using Jacobi Method Results
figure;
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,T_Jacobi,':','ShowText','on'),colorbar,colormap(jet)
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y')
title({['Temperature Contour Plot'];['Jacobi Method'];['Iterations taken : ',num2str(n_Jacobi)]})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. Gauss-Seidel Method
% Initializing Temperature Matrix
T_GS = zeros(N,M);

% Boundary Conditions
T_GS(N,:) = 900; % Bottom Boundary
T_GS(1,:) = 600; % Top Boundary
T_GS(:,1) = 400; % Left Boundary
T_GS(:,M) = 800; % Right Boundary

% Computation - Gauss Seidel Method
epsilon = 1e-8; % Tolerance or Error Sensitivity
err_GS = 1; % Error
Err_GS = []; % Error Array
n_GS = 0; % Iterator

while err_GS>epsilon
    err_GS = 0;
    for j=2:N-1
        for i=2:M-1
            T_GS_temp = T_GS(j,i);
            T_GS(j,i) = (1/ap)*(ae*T_GS(j+1,i) + aw*T_GS(j-1,i) + an*T_GS(j,i+1) + as*T_GS(j,i-1));
            err_GS = err_GS + power(T_GS(j,i)-T_GS_temp,2);
        end
    end
    err_GS = sqrt(err_GS/(M*N));
    Err_GS(end+1) = err_GS;
    n_GS = n_GS+1; % Iterator
end

% Plotting using Gauss-Seidel Method Results
figure;
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,T_GS,':','ShowText','on'),colorbar,colormap(jet)
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y')
title({['Temperature Contour Plot'];['Gauss-Seidel Method'];['Iterations taken : ',num2str(n_GS)]})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. PSOR Method
% Initializing Temperature Matrix
T_PSOR = zeros(N,M);
T_old_PSOR = zeros(N,M);

% Boundary Conditions
T_PSOR(N,:) = 900; % Bottom Boundary
T_PSOR(1,:) = 600; % Top Boundary
T_PSOR(:,1) = 400; % Left Boundary
T_PSOR(:,M) = 800; % Right Boundary

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
    for j=1:N
        for i=1:M
            T_old_PSOR(j,i) = T_PSOR(j,i);
        end
    end
    for j=2:N-1 % For Internal Grid points
        for i=2:M-1 % For Internal Grid points
            % T(j,i)_(k+1)th_iteration Calculation
            T_PSOR(j,i) = (1-w)*T_old_PSOR(j,i)+(w/(2*(1+B^2)))*((B^2)*T_PSOR(j,i-1)+T_PSOR(j-1,i)+T_old_PSOR(j,i+1)+(B^2)*T_old_PSOR(j+1,i));
        end
    end
    err_PSOR = 0; % Error
    for j=1:N
        for i=1:M
            err_PSOR = err_PSOR + power(T_PSOR(j,i)-T_old_PSOR(j,i),2);
        end
    end
    err_PSOR = sqrt(err_PSOR/(M*N));
    Err_PSOR(end+1) = err_PSOR;
    n_PSOR = n_PSOR+1; % Iterator
end

% Plotting using PSOR Method Results
figure;
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,T_PSOR,':','ShowText','on'),colorbar,colormap(jet)
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y')
title({['Temperature Contour Plot'];['PSOR Method'];['Iterations taken : ',num2str(n_PSOR)]})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Errors VS Iterations Comparison Plot
figure;
x_Jacobi = 1:1:n_Jacobi;
x_GS = 1:1:n_GS;
x_PSOR = 1:1:n_PSOR;
loglog(x_Jacobi,Err_Jacobi,'r-',x_GS,Err_GS,'b-',x_PSOR,Err_PSOR,'g-')
hold on;
x = 1:1:max([n_Jacobi,n_GS,n_PSOR]);
Tolerance = (1e-8)*ones(1,max([n_Jacobi,n_GS,n_PSOR]));
plot(x,Tolerance,'k-');
hold off;
xlabel('log(Iterations)'),ylabel('log(Error)')
title('Error VS Iterations Comparison')
legend('Jacobi','Gauss-Seidel','PSOR','Tolerance')
