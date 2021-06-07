% 2D Unsteady Heat Conduction Example-1
clear all; clc

% Known Parameters
Lx = 1; % Length of x-domain (1 m)
Ly = 1; % Length of y-domain (1 m)
alpha = 1.4; % Thermal Diffusivity
M = 41; % Number of grid points along x-axis
N = 41; % Number of grid points along y-axis

% Calculated Constants
dx = Lx/(M-1); % Grid Size in x-direction
dy = Ly/(N-1); % Grid Size in y-direction
dt = 1e-3; % Time Step
gamma_x = alpha*dx/(dt^2);
gamma_y = alpha*dy/(dt^2);

% Initializing Temperature matrix
T = zeros(N,M);
T_array = {}; % Temperature Array to store T at each time step

% Boundary Conditions
T(:,1) = 400; % Left Wall Boundary Condition
T(N,:) = 900; % Bottom Wall Boundary Condition
T(:,M) = 800; % Right Wall Boundary Condition
T(1,:) = 600; % Top Wall Boundary Condition

% Computation
n = 0; % Iterator
t = 0;
while t<=10
    err = 0; % Error
    for j=2:N-1
        for i=2:M-1
            T_temp = T(j,i);
            T(j,i) = (1/(1+2*gamma_x+2*gamma_y))*(gamma_x*(T(j+1,i)+T(j-1,i))+gamma_y*(T(j,i+1)+T(j,i-1))+T_temp);
            err = err + power(T(j,i)-T_temp,2);
        end
    end
    err = sqrt(err/(M*N));
    %T_array = [T_array;T];
    n = n+1;
    if n==10 || n==30 || n==50 || n==100 || n==150 || n==200
       T_array = [T_array;T];
    elseif n==275 || n==350 || n==550 || n==750 || n==1000 || n==1250
       T_array = [T_array;T]; 
    end
    t = t+dt;
end
%T
n
err

% Plotting
figure;
x = linspace(0,Lx,M);
y = linspace(Ly,0,N);
[X,Y] = meshgrid(x,y);
contourf(X,Y,T,':','ShowText','on'),colorbar,colormap(jet)
set(gca, 'XTick',0:1:Lx)
xlabel('X'),ylabel('Y'),title('Temperature Variation (T)')

% Animation
figure;
myVideo = VideoWriter('Unsteady Heat COnduction'); %open video file
myVideo.FrameRate = 5;
open(myVideo)

for k=1:12
    %figure;
    x = linspace(0,Lx,M);
    y = linspace(Ly,0,N);
    [X,Y] = meshgrid(x,y);
    TT = T_array(k);
    contourf(X,Y,TT{:,:},':','ShowText','on'),colorbar,colormap(jet)
    xlabel('X'),ylabel('Y'),title('Temperature Variation (T)')
    pause(0.02)
    frame = getframe(gcf);
    writeVideo(myVideo,frame);
end
close(myVideo)
