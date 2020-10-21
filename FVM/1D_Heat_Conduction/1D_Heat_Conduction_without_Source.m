%Steady State 1D Heat Conduction without Source -> using Finite Volume Method (FVM)

clear all; clc;

%Following SI Units
%Known Values
L = 1; %Length of the domain (or rod in this example)
n = 5; %Number of control volumes considered
dx = L/n; %Grid spacing; %Assumed Uniform Grid Spacing

%Thermal Conductivity
%Assumed constant throught the domain
ke = 1000;
kw = 1000;

%Boundary Conditions
%Temperature in Degree Celsius
TA = 100;
TB = 500;

%Calculated Constants
aw = kw/dx;
ae = ke/dx;

%Mesh Creation
x = dx/2:dx:1;
N = length(x);

%Matrices creation
A = zeros(N,N); %Coefficient Matrix
B = zeros(N,1); %Constants' Matrix

%Algorithm
for i=1:N
    if x(i) == min(x) %For Boundary Node-1 (Left most Node)
        A(i,i) = ae + 2*(kw/dx);
        A(i,i+1) = -ae;
        b(i,1) = 2*(kw/dx)*TA;
    elseif x(i) == max(x) %For Boundary Node-2 (Right most Node)
        A(i,i) = aw + 2*(ke/dx);
        A(i,i-1) = -aw;
        b(i,1) = 2*(ke/dx)*TB;
    else %For all Internal Nodes
        A(i,i) = ae + aw;
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
    end
end

%Final results
A %Coefficient Tri-Diagonal Matrix
b %Constant Matrix
T = A\b %Displays the Temperature matrix after solving

%Plots and Comparison between Analytical and Exact solution curves
y = linspace(0,1,500);
Te = ((TB-TA)/L)*y + TA; %Exact solution for Temperature
plot(x,T,'o-',y,Te,'r')
xlabel('Length')
ylabel('Temperature')
legend('Analytical solution','Exact solution')
