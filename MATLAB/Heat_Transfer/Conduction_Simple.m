clear all
close all
clc

Lx = 1.0;
Ly = 1.0;

Nx = 100;
Ny = 100;

dx = Lx/Nx;
dy = Ly/Ny;

x = linspace(0,Lx,Nx+1);
y = linspace(0,Ly,Ny+1);

[X,Y]=meshgrid(x,y);
X=X';
Y=Y';

T = zeros(Nx+2,Ny+2);
T = T+(20+273);

t = 300.0;

alpha = 10^-4;
k = 205.0;
g = 0;

dt = (1/(2*alpha))*(1/((1/(dx^2))+(1/(dy^2))));

nt = ceil(t/dt);

tstep = 0;

% Boundary Condition Parameters

q1 = 0;
q2 = 0;

T1 = 100+273;
T2 = 100+273;

while (tstep < nt)
    
    T_old = T;
    

    T(2:Nx+1,2:Ny+1)= T_old(2:Nx+1,2:Ny+1)+(((T_old(2:Nx+1,3:Ny+2)+T_old(2:Nx+1,1:Ny)-T_old(2:Nx+1,2:Ny+1).*2)/(dy*dy))...\
                                           +((T_old(3:Nx+2,2:Ny+1)+T_old(1:Nx,2:Ny+1)-T_old(2:Nx+1,2:Ny+1).*2)/(dx*dx)))*alpha*dt + (dt/k)*g;
    
      
                                       
    T(1,:) = T(2,:) + (q1*dx)/k;
    T(:,1) = T(:,2) + (q2*dx)/k;

    T(Nx+2,:) = 2*T1 - T(Nx+1,:);
    T(:,Ny+2) = 2*T2 - T(:,Ny+1);
    
 
    T_mid = T(Nx/2,Ny/2);
    
    if(T_mid >= 70+273)
        exit
    end
    %T_res = sum(sum((T-T_old).^2))/((Nx+2)*(Ny+2));
    
    %if (T_res < .000001 && T_res ~=0) 
    %    exit
    %end
    
    tstep = tstep + 1;
    
end

TT = (T(1:Nx+1,1:Ny+1)+T(1:Nx+1,2:Ny+2)+T(2:Nx+2,1:Ny+1)+T(2:Nx+2,2:Ny+2))/4;

figure
contourf(X,Y,TT-273,'LineStyle','none')
figure
contour(X,Y,TT-273,[linspace(0,50,50),linspace(50,350,50)])
