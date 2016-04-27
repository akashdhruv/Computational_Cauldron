% Matlab Script to Simulate Conduction over a 2D - Flate Plate.
% MAE 3187 - Heat Transfer
% Author - Jon Snow

% Emptying the stack memory from previous scripts, closing all figures and
% clearing screen
clear all
close all
clc

% All values are in SI Units

% Defining Domain, Lx = 30 cm, Ly = 30 cm
Lx = 0.3; % m
Ly = 0.3; % m

% Number of grid points
Nx = 100;
Ny = 100;

% Grid spacing
dx = Lx/Nx;
dy = Ly/Ny;

% Defining node points
x = linspace(0,Lx,Nx+1);
y = linspace(0,Ly,Ny+1);

% Creating grid - will be useful when plotting - See post
% processing in the end
[X,Y]=meshgrid(x,y);
X=X';
Y=Y';

% Defining Temperature points on cell centers not on nodes points. You will
% see why.
T = zeros(Nx+2,Ny+2);

% Setting initial values at room temperature
T = T+(40+273);

% Simulation time
t = 90.0;

% Thermal Conductivity
k = 205; % W/mK

% Thermal Diffusivity
alpha = 8.344*(10^-5); % m^2/s

% Heat Generation
g = 10^6; % W/m^3

% Calculating dt using stability condition
dt = 0.5*(1/(2*alpha))*(1/((1/(dx^2))+(1/(dy^2))));

% Defining number of time steps
nt = ceil(t/dt); 

%nt = floort(t/dt); % can also use this

% ceil (ceiling) and floor (floor) return high/low integer values
% respectively. Type ' help ceil ' or ' help floor ' in command window 
% to learn more about them

tstep = 0; % time loop counter. Start with 0

% Boundary Condition Parameters

q1 = 10^6; % Flux - Low x boundary condition
q2 = 5*(10^5); % Flux - Low y boundary condition

T1 = 10+273; % Temp - High x boundary condition
T2 = 10+273; % Temp - High y boundary condition

% Begin Time Loop
while (tstep < nt)
    
    T_old = T; % Store T values from previous time step in T_old
    
    
    % Using vectorized equation to apply governing equation at INNER NODES
    % Eliminates the use of two for loops - increases computational speed
    
    T(2:Nx+1,2:Ny+1)= T_old(2:Nx+1,2:Ny+1)+...\
                   (((T_old(2:Nx+1,3:Ny+2)+T_old(2:Nx+1,1:Ny)-T_old(2:Nx+1,2:Ny+1).*2)/(dy*dy))...\
                   +((T_old(3:Nx+2,2:Ny+1)+T_old(1:Nx,2:Ny+1)-T_old(2:Nx+1,2:Ny+1).*2)/(dx*dx)))*alpha*dt...\
                    +(alpha*dt/k)*g;
     
    % Applying BCs. Using cell center method makes BC application simple.
    % These points are not on nodes but outside computational domain. These
    % are called Ghost Cells or Guard Cells.
    
    % Low x BC
    T(1,:) = T(2,:) + (q1*dx)/k;
    
    % Low y BC
    T(:,1) = T(:,2) + (q2*dy)/k;

    % High x BC
    T(Nx+2,:) = 2*T1 - T(Nx+1,:);
    
    % High y BC
    T(:,Ny+2) = 2*T2 - T(:,Ny+1);

    % Calculating residuals 
    % Error = RMS(T - T_old)
    
    T_res = sum(sum((T-T_old).^2))/((Nx+2)*(Ny+2));
    
    % If Error < Tolerance , then exit the time loop
    
    if (T_res < .000001 && T_res ~=0) 
        break
    end
    
    % increment time counter
    tstep = tstep + 1;
    
end
% End Time Loop

% Post Processing

% Interpolate values from cell center to nodal points
TT = (T(1:Nx+1,1:Ny+1)+T(1:Nx+1,2:Ny+2)+T(2:Nx+2,1:Ny+1)+T(2:Nx+2,2:Ny+2))/4;

% Plot results, use X,Y meshgrid points to specify T location
figure
contourf(X,Y,TT-273,'LineStyle','none')