%Point Gauss Seidel Algorithm for 2D Steady State Heat Conduction

clear all; clc;
nx = 6; %Number of nodes along x-direction
ny = 6; %Number of nodes along y-direction

%Length of the domain (SI Units)
Lx = 1;
Ly = 1;

%Grid Spacing
dx = 0.2; %Delta_X
dy = 0.2; %Delta_Y
B = dx/dy; %Beta = (Delta_X/Delta_Y)

%Initializing Temperature Matrix with all values as zeros
T = zeros(nx,ny);

%Boundary Conditions
T(1,:) = 100;
T(:,1) = 100;

%Algorithm
error = 1; %Error Aprroximation; %Can be arbitrary
epsilon = 0.001; %Allowable Tolerance
N = 0; %Iteration Counter

while error>=epsilon %This step goes on as long as Max error in the "T matrix" is greater than or equal to the Tolerance (Epsilon)

  T_old = T; %T_old refers to the T obtained in the previous iteration.
  
  N = N+1; %Iteration counter incremented for every iteration that the "while loop" takes
  
  for i=2:nx-1 %This double "for loop" runs for obtaining all the Temperature values for all the "INTERNAL NODES".
    for j=2:ny-1
      T(i,j) = (1/(2*(1+beta^2)))*(T_old(i+1,j) + T(i-1,j) + (beta^2)*(T_old(i,j+1) + T(i,j-1)));
    end
  end
  
  error = max(max(abs(T-T_old)));
end

%Results
T %Displays the resulting Temperature Matrix
N % No. of iterations taken
error %Final error

%Contour plotting
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);
contour(X,Y,T)
contour(X,Y,T,'Fill','on'),colorbar %can also use contourf(X,Y,T)
xlabel('X'),ylabel('Y'),title('Temperature variation')
