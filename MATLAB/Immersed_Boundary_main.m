% Code to compute the Navier-Stokes Solution for a Lid-Driven Cavity

clear all
close all
clc

Ly=1; % Domain Lengths
Lx=1;


mu=.01; % Viscosity for Re=100

U=1; % Velocity
w=1.1; % SOR Relaxation Factor

Nx=50; % Grid Size
Ny=50;

dx=Lx/(Nx-1); % dx and dy
dy=Ly/(Ny-1);

% Defining Staggered grid

xstart=0;
ystart=0;

x=zeros(Nx,1);
y=zeros(Ny,1);

for i=1:Nx
x(i)=xstart+(i-1)*dx;
end

for i=1:Ny
y(i)=ystart+(i-1)*dy;
end

% Pressure points
xp=zeros(Nx+1,1);
yp=zeros(Ny+1,1);
p=zeros(Nx+1,Ny+1);

% u-velocity points
xu=zeros(Nx,1);
yu=zeros(Ny+1,1);
u=zeros(Nx,Ny+1);
ut=zeros(Nx,Ny+1); % Predicted Velocity Array
G1=zeros(Nx-1,Ny); % Will be used in RK3

% v-velocity points
xv=zeros(Nx+1,1);
yv=zeros(Ny,1);
v=zeros(Nx+1,Ny);
vt=zeros(Nx+1,Ny); % Predicted Velocity Array
G2=zeros(Nx,Ny-1); % Will be used in RK3


% Defining Pressure Points
for i=2:Nx
    for j=2:Ny
        xp(i)=(x(i-1)+x(i))/2;
        yp(j)=(y(j-1)+y(j))/2;
    end
end
xp(1)=xp(2)-dx;
yp(1)=yp(2)-dy;

xp(end)=xp(end-1)+dx;
yp(end)=yp(end-1)+dy;

% Defining u-velocity points
for i=1:Nx
    for j=2:Ny
        xu(i)=x(i);
        yu(j)=(y(j-1)+y(j))/2;
    end
end
xu(1)=xu(2)-dx;
yu(1)=yu(2)-dy;

xu(end)=xu(end-1)+dx;
yu(end)=yu(end-1)+dy;

% Defining v-velocity points

for i=2:Nx
    for j=1:Ny
        xv(i)=(x(i-1)+x(i))/2;
        yv(j)=y(j);
    end
end

xv(1)=xv(2)-dx;
yv(1)=yv(2)-dy;

xv(end)=xv(end-1)+dx;
yv(end)=yv(end-1)+dy;


% Plotting the Grid
figure
hold on
mesh(x,y,zeros(Nx,Ny)) % Mesh

% for i=1:Nx+1
%     for j=1:Ny+1
%         plot(xp(i),yp(j),'.b')
%     end
% end
% for i=1:Nx
%     for j=1:1
%         plot(x(i),y(j),'.r')
%         plot(x(i),y(Ny),'.r')
%     end
% end
% for i=1:Nx
%     for j=1:Ny+1
%         plot(xu(i),yu(j),'xk')
%     end
% end
% for i=1:Nx+1
%     for j=1:Ny
%         plot(xv(i),yv(j),'*m')
%     end
% end

for i=20:30
    xb1(i-19)=x(i);
    yb1(i-19)=y(i);
end
xb=[fliplr(xb1),xb1(1)*ones(1,10),xb1,xb1(end)*ones(1,10)];
yb=[yb1(1)*ones(1,10),yb1,yb1(end)*ones(1,10),fliplr(yb1)];
plot(xb,yb,'k')
fill(xb,yb,'k')
%% Unsteady N-S Solver

t=0; % t-initial
tmax=10; % t-final
nt=900; %  number of time steps, determined using CFL condition
dt=(tmax-t)/nt; % dt

[X Y]=meshgrid(x,y);

figure

% Starting time loop
for tstep=1:nt
    
    % RK3 First Step
    
    % Applying Velocity Boundary Conditions
    u(:,1)=-u(:,2);
    u(:,end)=2*U-u(:,end-1);
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);

    % Applying Pressure Boundary Conditions
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);
    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    
    % Calculating Predicted Velocities
    for i=2:Nx-1
        for j=2:Ny
            if (i>=20 && i<=30 && j>=20 && j<=30)
                ut(i,j)=0;
            else
            ue=(u(i,j)+u(i+1,j))/2;
            uw=(u(i,j)+u(i-1,j))/2;
            us=(u(i,j)+u(i,j-1))/2;
            un=(u(i,j)+u(i,j+1))/2;
            vs=(v(i,j-1)+v(i+1,j-1))/2;
            vn=(v(i,j)+v(i+1,j))/2;
        
            G1(i,j)=(F1C(ue,uw,us,un,vs,vn,dx,dy)+FV(u(i,j),u(i+1,j),u(i-1,j),u(i,j+1),u(i,j-1),dx,dy,mu));
            ut(i,j)=u(i,j)+(dt/3)*G1(i,j); 
            end
        end
    end

    for i=2:Nx
        for j=2:Ny-1
            if (i>=20 && i<=30 && j>=20 && j<=30)
                vt(i,j)=0;
            else
        vn=(v(i,j)+v(i,j+1))/2;
        vs=(v(i,j)+v(i,j-1))/2;
        ve=(v(i,j)+v(i+1,j))/2;
        vw=(v(i,j)+v(i-1,j))/2;
        ue=(u(i,j)+u(i,j+1))/2;
        uw=(u(i-1,j)+u(i-1,j+1))/2;
        
        G2(i,j)=(F2C(vn,vs,ve,vw,ue,uw,dx,dy)+FV(v(i,j),v(i+1,j),v(i-1,j),v(i,j+1),v(i,j-1),dx,dy,mu));
        vt(i,j)=v(i,j)+(dt/3)*G2(i,j);
            end
        end
    end
    
    % Boundary Conditions for Predicted Velocities
    ut(:,1)=-ut(:,2);
    ut(:,end)=2*U-ut(:,end-1);
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);
    
    

    % Solving Poisson Equation 
    maxiter=500;

    for it=1:maxiter
        for i=2:Nx
            for j=2:Ny
                p(i,j)=w*0.25*(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-(3*dx/dt)*(vt(i,j)-vt(i,j-1)+ut(i,j)-ut(i-1,j)))+(1-w)*p(i,j);    
            end
        end
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);

    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    end


    % Velocities after first RK-3 Step
    for i=2:Nx-1
        for j=2:Ny
        u(i,j)=ut(i,j)-(dt/(dx*3))*(p(i+1,j)-p(i,j));
        end
    end
    for i=2:Nx
        for j=2:Ny-1
        v(i,j)=vt(i,j)-(dt/(dx*3))*(p(i,j+1)-p(i,j));
        end
    end
    
    %RK-3 Second Step
    
    % Applying Boundary Conditions for Velocites
    u(:,1)=-u(:,2);
    u(:,end)=2*U-u(:,end-1);
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);

    % Calculating Predicted Velocities
    
    for i=2:Nx-1
        for j=2:Ny
            if (i>=20 && i<=30 && j>=20 && j<=30)
                ut(i,j)=0;
            else
            ue=(u(i,j)+u(i+1,j))/2;
            uw=(u(i,j)+u(i-1,j))/2;
            us=(u(i,j)+u(i,j-1))/2;
            un=(u(i,j)+u(i,j+1))/2;
            vs=(v(i,j-1)+v(i+1,j-1))/2;
            vn=(v(i,j)+v(i+1,j))/2;
        
            G1(i,j)=-(5/9)*G1(i,j)+(F1C(ue,uw,us,un,vs,vn,dx,dy)+FV(u(i,j),u(i+1,j),u(i-1,j),u(i,j+1),u(i,j-1),dx,dy,mu));
            ut(i,j)=u(i,j)+(15*dt/16)*G1(i,j); 
            end
        end
    end

    for i=2:Nx
        for j=2:Ny-1
            if (i>=20 && i<=30 && j>=20 && j<=30)
                vt(i,j)=0;
            else
            vn=(v(i,j)+v(i,j+1))/2;
            vs=(v(i,j)+v(i,j-1))/2;
            ve=(v(i,j)+v(i+1,j))/2;
            vw=(v(i,j)+v(i-1,j))/2;
            ue=(u(i,j)+u(i,j+1))/2;
            uw=(u(i-1,j)+u(i-1,j+1))/2;
        
            G2(i,j)=-(5/9)*G2(i,j)+(F2C(vn,vs,ve,vw,ue,uw,dx,dy)+FV(v(i,j),v(i+1,j),v(i-1,j),v(i,j+1),v(i,j-1),dx,dy,mu));
            vt(i,j)=v(i,j)+(15*dt/16)*G2(i,j);
            end
        end
    end

    % Boundary Condition for Predicted Velocities
    ut(:,1)=-ut(:,2);
    ut(:,end)=2*U-ut(:,end-1);
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);

    % Solving Poisson Equation
    for it=1:maxiter
        for i=2:Nx
            for j=2:Ny
                p(i,j)=w*0.25*(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-(12*dx/(5*dt))*(vt(i,j)-vt(i,j-1)+ut(i,j)-ut(i-1,j)))+(1-w)*p(i,j);    
            end
        end
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);

    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    end

    % Velocities after Second RK-3 Step
    for i=2:Nx-1
        for j=2:Ny
            u(i,j)=ut(i,j)-(5*dt/(dx*12))*(p(i+1,j)-p(i,j));
        end
    end
    for i=2:Nx
        for j=2:Ny-1
            v(i,j)=vt(i,j)-(5*dt/(dx*12))*(p(i,j+1)-p(i,j));
        end
    end
    
    % RK-3 Third Step
    
    % Applying Boundary Conditions for Velocities
    u(:,1)=-u(:,2);
    u(:,end)=2*U-u(:,end-1);
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);


    % Calculating Predicted Velocities
    for i=2:Nx-1
        for j=2:Ny
            if (i>=20 && i<=30 && j>=20 && j<=30)
                ut(i,j)=0;
            else
        ue=(u(i,j)+u(i+1,j))/2;
        uw=(u(i,j)+u(i-1,j))/2;
        us=(u(i,j)+u(i,j-1))/2;
        un=(u(i,j)+u(i,j+1))/2;
        vs=(v(i,j-1)+v(i+1,j-1))/2;
        vn=(v(i,j)+v(i+1,j))/2;
        
        G1(i,j)=-(153/128)*G1(i,j)+(F1C(ue,uw,us,un,vs,vn,dx,dy)+FV(u(i,j),u(i+1,j),u(i-1,j),u(i,j+1),u(i,j-1),dx,dy,mu));
        ut(i,j)=u(i,j)+(8*dt/15)*G1(i,j);    
            end
        end
    end

    for i=2:Nx
        for j=2:Ny-1
            if (i>=20 && i<=30 && j>=20 && j<=30)
                vt(i,j)=0;
            else
            vn=(v(i,j)+v(i,j+1))/2;
            vs=(v(i,j)+v(i,j-1))/2;
            ve=(v(i,j)+v(i+1,j))/2;
            vw=(v(i,j)+v(i-1,j))/2;
            ue=(u(i,j)+u(i,j+1))/2;
            uw=(u(i-1,j)+u(i-1,j+1))/2;
        
            G2(i,j)=-(153/128)*G2(i,j)+(F2C(vn,vs,ve,vw,ue,uw,dx,dy)+FV(v(i,j),v(i+1,j),v(i-1,j),v(i,j+1),v(i,j-1),dx,dy,mu));
            vt(i,j)=v(i,j)+(8*dt/15)*G2(i,j);
            end
        end
    end

    % Boundary Conditions for the Predicted Velocities
    ut(:,1)=-ut(:,2);
    ut(:,end)=2*U-ut(:,end-1);
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);


    % Solving Poisson Equation
    for it=1:maxiter
        for i=2:Nx
            for j=2:Ny
                p(i,j)=w*0.25*(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-(4*dx/dt)*(vt(i,j)-vt(i,j-1)+ut(i,j)-ut(i-1,j)))+(1-w)*p(i,j);    
            end
        end
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);

    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    end

    % Calculating the Velocities after the final RK-3 step
    for i=2:Nx-1
        for j=2:Ny
        u(i,j)=ut(i,j)-(dt/(dx*4))*(p(i+1,j)-p(i,j));
        end
    end
    for i=2:Nx
        for j=2:Ny-1
        v(i,j)=vt(i,j)-(dt/(dx*4))*(p(i,j+1)-p(i,j));
        end
    end
    
    %RK-3 Complete


% Calculating the Velocities and Vorticity at Grid points

uu(1:Nx,1:Ny)=0.5*(u(1:Nx,2:Ny+1)+u(1:Nx,1:Ny));
vv(1:Nx,1:Ny)=0.5*(v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny));
ww(1:Nx,1:Ny)=(u(1:Nx,2:Ny+1)-u(1:Nx,1:Nx)-v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny))/(2*dx);

% Plotting Velocity Vectors and Vorticity Contours
hold off,quiver(flipud(rot90(uu)),flipud(rot90(vv)),'k');hold on;
contour(flipud(rot90(ww)),100),axis equal; axis([1 Nx 1 Ny]);
fill(xb,yb,'w')
%contour(flipud(rot90(p)),100),axis equal; axis([1 Nx 1 Ny]);

drawnow
pause(0.01)

t=t+dt % Advancing Time Step
end

%% Post Processing

% Velocity Potential and Stream Function Calculation
[phi,psi]=flowfun(uu',vv');

% Plotting Streamfunction Contours
figure
contourf(X,Y,psi); hold on; axis equal

% Plotting Streamlines and Velocity Vectors
figure
hold on
plot([0,0,1.0,1.0,0],[0,1.0,1.0,0,0],'k');
quiver(X,Y,uu',vv',2,'b')
contour(X,Y,psi,'k')
axis([0 1.0 0 1.0])
fill(xb,yb,'w');
axis equal

