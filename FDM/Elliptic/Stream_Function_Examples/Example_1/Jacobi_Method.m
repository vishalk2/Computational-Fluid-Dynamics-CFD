% Stream Function Example-1 using Jacobi Method

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
err = 1; % Error
n = 0; % Iterator

while err>epsilon
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
    err = 0; % Error
    for j=1:N
        for i=1:M
            err = err + power(psi_new(j,i)-psi_old(j,i),2);
        end
    end
    err = sqrt(err/(M*N));
    n = n+1; % Iterator
end
psi_new
n
err

% Plotting
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,psi_new,'ShowText','on'),colorbar
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y'),title('Stream Function(PSI) Contour Plot')
