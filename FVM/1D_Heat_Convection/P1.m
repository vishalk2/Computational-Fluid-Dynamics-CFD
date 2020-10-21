% 1D Heat Convection & Diffusion -> using Finite Volume Method (FVM)

clear all; clc;

%Following SI Units

%Known Values
L = 1; %Length of the domain
n = 5; %Number of Control Volumes
dx = L/n; %Grid Spacing; %Assumed Uniform grid spacing
rho = 1; %Density of the bulk
gamma = 0.1; %Viscocity
u = 0.1; %Velocity of the bulk fluid

%Boundary Conditions
%Phi => some parameter
PhiA = 1; %At left boundary
PhiB = 0; %At right boundary

%Calculated constants
Dw = Gamma/dx; De = Gamma/dx;
Fw = rho * u; Fe = rho * u;
Area = dx * dx;
aw = Dw + (Fw/2); ae = De - (Fe/2);

%Create Mesh
x = dx/2:dx:L;
N = length(x);

%Matrix creation
A = zeros(N,N);
b = zeros(N,1);

%Algorithm
for i = 1:N
    if x(i) == min(x) %For Boundary Point-1 (left boundary)
        A(i,i) = ae + ((2*Dw)+Fw) + (Fe-Fw);
        A(i,i+1) = -ae;
        b(i,1) = ((2*Dw)+Fw)*PhiA;
    elseif x(i) == max(x) %For Boundary Point-2 (right boundary)
        A(i,i) = aw + ((2*De)-Fe) + (Fe-Fw);
        A(i,i-1) = -aw;
        b(i,1) = ((2*De)-Fe)*PhiB;
    else %For Internal Nodes
        A(i,i) = ae + aw + (Fe - Fw);
        A(i,i+1) = -ae;
        A(i,i-1) = -aw;
    end
end

%Final Results
%Final results
A %Coefficient Tri-Diagonal Matrix
b %Constant Matrix
Phi = A\b %Displays the Phi matrix after solving

%Plots and comparison between Analytical/Exact & Numerical solution curve
y=linspace(0,1,500);
Phi_e=1-((exp(rho*u*y/Gamma)-1)/(exp(rho*u*L/Gamma)-1)); %Exact solution (Analytical) for Phi
plot(x,Phi,'o-',y,Phi_e,'r')
xlabel('length')
ylabel('Phi')
legend('Numerical value','Exact/Analytical value')

%Below code gives only the plot of the curve obtained from Numerival solution
% plot([0,x,L],[PhiA;Phi;PhiB],'o-b')
% title('Property Transported by Convetion and Diffusion')
% xlabel('meters (m)')
% ylabel('Phi')
