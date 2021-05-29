% Temperature Variation Example-4 using PSOR Method

% Known Parameters
Lx = 1; % Length of domain along x-direction
Ly = 1; % Length of domain along y-direction
M = 41; % Number of grid points starting from 1 in x-direction
N = 41; % Number of grid points starting from 1 in y-direction

% Calculated constants
dx = Lx/(M-1); % Grid size along x-direction
dy = Ly/(N-1); % Grid size along y-direction

% Finding Optimum Relaxation Factor
B = dx/dy;
a = ((cos(pi/(M-1))+(B^2)*cos(pi/(N-1)))/(1+(B^2)))^2;
w = (2-sqrt(1-a))/a; % Optimum Relaxation Factor

% Initializing Temperature Matrix
T = zeros(N,M);

% Boundary Conditions
T(N,:) = 100; % Bottom Boundary T = 100 C
T(1,:) = 0; % Top Boundary T = 0 C
T(:,1) = 100; % Left Boundary T = 100 C
T(:,M) = 0; % Right Boundary T = 0 C

% Computation
epsilon = 1e-8; % Error Sensivity
err = 1; % Error
n = 0; % Iterator

while err>epsilon
    for j=1:N
        for i=1:M
            T_old(j,i) = T(j,i);
        end
    end
    for j=2:N-1 % For Internal Grid points
        for i=2:M-1 % For Internal Grid points
            % T(j,i)_(k+1)th_iteration Calculation
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
w
T
n
err

% Plotting
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,T,'ShowText','on'),colorbar
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y'),title('Temperature Contour Plot')
