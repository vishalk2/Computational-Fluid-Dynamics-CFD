%Gauss Seidel Algorithm & Point Successive Over-Relaxation Method (PSOR) for Steady State 1D Heat Conduction

n = 15; %Defining number of nodes
L = 1; %Length of the domain (SI Units)

N = 0; %Iteration Counter => gives the total number of iterations taken to arrive at the approximate answer

%Initializing the Temperature matrix
T = zeros(1,n); %Initializes a "1 x n" matrix for temperature

%Boundary Conditions
%TA = 100; => Left Boundary Node
%TB = 0;   => Right Boundary Node
T(1,1) = 100;
T(1,n) = 0;

%Relaxation Parameter (Used for PSOR)
w = 1.8; %In general 1<w<2
%If w = 1, we arrive at the solution based on Gauss Seidel Method but with a lot more number of iterations

%Algorithm
error = 1; %Error initially; %Can be arbitrary
epsilon = 0.000001; %Allowable Tolerance
%IF ERROR <= EPSILON, we arrive at our solution for "N" number of iterations at the end

while error>=epsilon %This step goes on as long as Max error in the "T matrix" is greater than or equal to the Tolerance (Epsilon)

  T_old = T; %T_old refers to the T obtained in the previous iteration.
  
  N = N+1; %Iteration counter incremented for every iteration that the "while loop" takes
    
  for i=2:n-1 %This loop runs for obtaining all the Temperature values for all the "INTERNAL NODES".
    T(i) = (1-w)*T_old(i) + (w/2)*(T(i-1) + T_old(i+1));
  end
  
  error = max(abs(T-T_old));
end

%Results
T %Displays the approximated Temperature Matrix
N %Displays the number of iterations took to arrive at the final values
  
%Plots
x = linspace(1,L,n);
plot(x,T,'LineWidth',2)
xlabel('Length')
ylabel('Temperature')
title('Temperature variation in a 1D Steady Heat Conduction Rod')
legend('T = f(x)')
