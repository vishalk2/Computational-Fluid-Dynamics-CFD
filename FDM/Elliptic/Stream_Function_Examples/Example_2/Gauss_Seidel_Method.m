% Stream Function Example-2 using Gauss-Seidel Method

% Known Parameters
Lx = 4; % Length of domain along x-direction
Ly = 4; % Length of domain along y-direction
M = 21; % Number of grid points starting from 1 in x-direction
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
psi = zeros(N,M);

% Boundary Conditions
psi(N,:) = 0; % Bottom Boundary
psi(1,:) = 1; % Top Boundary
psi(:,1) = 0; % Left Boundary
psi(:,M) = 0; % Right Boundary

% Computation
epsilon = 1e-8; % Error Sensivity
err = 1; % Error
n = 0; % Iterator

while err>epsilon
    err = 0; % Error
    for j=2:N-1 % For Internal Grid points
        for i=2:M-1 % For Internal Grid points
            psi_temp = psi(j,i);
            % Psi(j,i)_(k+1)th_iteration Calculation
            psi(j,i) = (1/ap)*(ae*psi(j+1,i) + aw*psi(j-1,i) + an*psi(j,i+1) + as*psi(j,i-1));
            err = err + power(psi(j,i)-psi_temp,2);
        end
    end
    err = sqrt(err/(M*N));
    n = n+1; % Iterator
end
psi
n
err

% Plotting
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,psi,'ShowText','on'),colorbar
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y'),title('Stream Function(PSI) Contour Plot')
