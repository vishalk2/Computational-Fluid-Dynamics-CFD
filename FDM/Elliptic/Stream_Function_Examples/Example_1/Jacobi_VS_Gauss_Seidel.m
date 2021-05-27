% Comparison of Error variation against Iteration w.r.t Jacobi
% and Gauss-Seidel Methods

% Known Parameters
Lx = 6; % Length of domain along x-direction
Ly = 4; % Length of domain along y-direction
M = 31; % Number of grid points starting from 1 in x-direction
N = 21; % Number of grid points starting from 1 in y-direction

% Calculated constants
dx = Lx/(M-1); % Grid size along x-direction
dy = Ly/(N-1); % Grid size along y-direction
% Coefficients
ap = 2*((1/(dx^2))+(1/(dy^2)));
ae = 1/(dx^2);
aw = 1/(dx^2);
an = 1/(dy^2);
as = 1/(dy^2);

% Jacobi Method
% Initializing Stream Function Matrix
psi_new = zeros(N,M);
psi_old = zeros(N,M);

% Boundary Conditions
psi_new(N,1:6) = 0; % Bottom Boundary
psi_new(N,7:M) = 100; % Bottom Boundary
psi_new(:,1) = 0; % Left Boundary
psi_new(1,:) = 0; % Top Boundary

% Computation
epsilon = 1e-8; % Error Sensivity
err_1 = 1; % Error
Err_1 = []; % Error Array
n_1 = 0; % Iterator

while err_1>epsilon
    for j=1:N
        for i=1:M
            psi_old(j,i) = psi_new(j,i);
        end
    end
    for j=2:N-1 % For Internal Grid points
        for i=2:M-1 % For Internal Grid points
            % Psi(j,i)_(k+1)th_iteration Calculation
            psi_new(j,i) = (1/ap)*(ae*psi_old(j+1,i) + aw*psi_old(j-1,i) + an*psi_old(j,i+1) + as*psi_old(j,i-1));
        end
    end
    for j=1:N % Right Neumann Boundary Condition
        psi_new(j,M) = psi_new(j,M-1);
    end
    err_1 = 0; % Error
    for j=1:N
        for i=1:M
            err_1 = err_1 + power(psi_new(j,i)-psi_old(j,i),2);
        end
    end
    err_1 = sqrt(err_1/(M*N));
    Err_1(end+1) = err_1;
    n_1 = n_1+1; % Iterator
end

% Gauss-Seidel Method
% Initializing Stream Function Matrix
psi = zeros(N,M);

% Boundary Conditions
psi(N,1:6) = 0; % Bottom Boundary
psi(N,7:M) = 100; % Bottom Boundary
psi(:,1) = 0; % Left Boundary
psi(1,:) = 0; % Top Boundary

err_2 = 1; % Error
Err_2 = []; % Error Array
n_2 = 0; % Iterator

while err_2>epsilon
    err_2 = 0; % Error
    for j=2:N-1 % For Internal Grid points
        for i=2:M-1 % For Internal Grid points
            psi_temp = psi(j,i);
            % Psi(j,i)_(k+1)th_iteration Calculation
            psi(j,i) = (1/ap)*(ae*psi(j+1,i) + aw*psi(j-1,i) + an*psi(j,i+1) + as*psi(j,i-1));
            err_2 = err_2 + power(psi(j,i)-psi_temp,2);
        end
    end
    for j=1:N % Right Neumann Boundary Condition
        psi(j,M) = psi(j,M-1);
    end
    err_2 = sqrt(err_2/(M*N));
    Err_2(end+1) = err_2;
    n_2 = n_2+1; % Iterator
end

% Outputs
display('Number of Iterations in Jacobi Method:',num2str(n_1));
display('Number of Iterations in Gauss-Seidel Method:',num2str(n_2));
%Err_1
%Err_2

% Comparison Plots
x1 = 1:1:n_1;
x2 = 1:1:n_2;
loglog(x1,Err_1,'r-',x2,Err_2,'b-')
xlabel('log(Iterations)'),ylabel('log(Error)')
title('Error VS Iterations Comparison')
legend('Jacobi','Gauss-Seidel')
