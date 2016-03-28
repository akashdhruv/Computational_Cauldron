% Code to compute the Navier-Stokes Solution for a Lid-Driven Cavity
% Written by Akash Dhruv
clear all
close all
clc

Ly=1; % Domain Lengths
Lx=1;


mu=.001; % Viscosity for Re=100

U=1; % Velocity
w=1.0; % SOR Relaxation Factor

Nx=30; % Grid Size
Ny=30;

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
%     for j=1:Ny
%         plot(x(i),y(j),'.r')
%     end
% end
% for i=1:Nx
%     for j=1:Ny+1
%         plot(xu(i),yu(j),'xc')
%     end
% end
% for i=1:Nx+1
%     for j=1:Ny
%         plot(xv(i),yv(j),'*m')
%     end
% end

% xcircle=linspace(0.8,.2,50);
% ycircle=sqrt(0.3^2-(xcircle-0.5).^2);
% 
% xcircle=[xcircle,fliplr(xcircle(1:end-1))];
% ycircle=[ycircle,-fliplr(ycircle(1:end-1))];
% 
% ycircle=ycircle+0.5;
% 
% plot(xcircle,ycircle,'k')
% 
% axis equal
% 
% egu=1;
% for i=1:length(xu)
%     for j=1:length(yu)
%         
%         if(sqrt((xu(i)-0.5)^2+(yu(j)-0.5)^2)<=0.3)
%             iu(egu,:)=[i,j];
%             egu=egu+1;
%         end
%         
%     end
% end
% for i=1:egu-1
%     plot(xu(iu(i,1)),yu(iu(i,2)),'.b')
% end
% egv=1;
% for i=1:length(xv)
%     for j=1:length(yv)
%         
%         if(sqrt((xv(i)-0.5)^2+(yv(j)-0.5)^2)<=0.3)
%             iv(egv,:)=[i,j];
%             egv=egv+1;
%         end
%         
%     end
% end
% for i=1:egv-1
%     plot(xv(iv(i,1)),yv(iv(i,2)),'xk')
% end
%xb=linspace(0,0.5,20);
%yb=0.2*ones(size(xb));
%yb=[yb,linspace(0.2,0,20)];
%xb=[xb,0.5*ones(1,20)];

%xbo=linspace(0.5,0,20);
%ybo=0*ones(size(xbo));
%ybo=[ybo,linspace(0,0.2,20)];
%xbo=[xbo,0*ones(1,20)];

%plot(xb,yb,'k')
%axis equal

%iv=find(xv<0.5);
%iv=iv(end);

%iu=find(xu<0.5);
%iu=iu(end);

%jv=find(yv<0.2);
%jv=jv(end);

%ju=find(yu<0.2);

%ju=ju(end);

%for i=1:iu
%    for j=ju:ju
%        plot(xu(i),yu(j),'.b')
%        plot(xu(i),yu(j+1),'xb')
%    end
%end

%for i=iu:iu
%    for j=2:ju
%        plot(xu(i),yu(j),'.r')
%        plot(xu(i+1),yu(j),'xr')
%    end
%end


%for i=2:iv
%    for j=jv:jv
%        plot(xv(i),yv(j),'.k')
%        plot(xv(i),yv(j+1),'xk')
%    end
%end

%for i=iv:iv
%    for j=1:jv
%        plot(xv(i),yv(j),'.m')
%        plot(xv(i+1),yv(j),'xm')
%    end
%end
axis equal
%%
% Unsteady N-S Solver

t=0; % t-initial
tmax=5; % t-final
dt=0.34*min(dx,dy);
nt=ceil((tmax-t)/dt);
maxiter=100;
[X Y]=meshgrid(x,y);

figure

% Starting time loop
for tstep=1:nt
    
    % RK3 First Step
    
     % Applying Boundary Conditions for Velocites
    u(:,1)=-u(:,2);
    u(:,end)=2-u(:,end-1);
    u(1,:)=0;
    u(end,:)=0;
 %   u(1:iu,ju)=-u(1:iu,ju+1);
 %   u(iu,2:ju)=-u(iu+1,2:ju);
    
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);
    v(:,1)=0;
    v(:,end)=0;
 %   v(iv,1:jv)=-v(iv+1,1:jv);
 %   v(2:iv,jv)=-v(2:iv,jv+1);

    % Calculating Predicted Velocities

    ue1=(u(2:Nx-1,2:Ny)+u(3:Nx,2:Ny))/2;
    uw1=(u(2:Nx-1,2:Ny)+u(1:Nx-2,2:Ny))/2;
    us1=(u(2:Nx-1,2:Ny)+u(2:Nx-1,1:Ny-1))/2;
    un1=(u(2:Nx-1,2:Ny)+u(2:Nx-1,3:Ny+1))/2;
    vs1=(v(2:Nx-1,1:Ny-1)+v(3:Nx,1:Ny-1))/2;
    vn1=(v(2:Nx-1,2:Ny)+v(3:Nx,2:Ny))/2;
        
    G1(2:Nx-1,2:Ny)=(F1C(ue1,uw1,us1,un1,vs1,vn1,dx,dy)+FV(u(2:Nx-1,2:Ny),u(3:Nx,2:Ny),u(1:Nx-2,2:Ny)...
                             ,u(2:Nx-1,3:Ny+1),u(2:Nx-1,1:Ny-1),dx,dy,mu));
    ut(2:Nx-1,2:Ny)=u(2:Nx-1,2:Ny)+(dt/3)*G1(2:Nx-1,2:Ny); 
    
     
    vn12=(v(2:Nx,2:Ny-1)+v(2:Nx,3:Ny))/2;
    vs12=(v(2:Nx,2:Ny-1)+v(2:Nx,1:Ny-2))/2;
    ve12=(v(2:Nx,2:Ny-1)+v(3:Nx+1,2:Ny-1))/2;
    vw12=(v(2:Nx,2:Ny-1)+v(1:Nx-1,2:Ny-1))/2;
    ue12=(u(2:Nx,2:Ny-1)+u(2:Nx,3:Ny))/2;
    uw12=(u(1:Nx-1,2:Ny-1)+u(1:Nx-1,3:Ny))/2;
        
    G2(2:Nx,2:Ny-1)=(F2C(vn12,vs12,ve12,vw12,ue12,uw12,dx,dy)+FV(v(2:Nx,2:Ny-1),v(3:Nx+1,2:Ny-1),...
                             v(1:Nx-1,2:Ny-1),v(2:Nx,3:Ny),v(2:Nx,1:Ny-2),dx,dy,mu));
    vt(2:Nx,2:Ny-1)=v(2:Nx,2:Ny-1)+(dt/3)*G2(2:Nx,2:Ny-1);
    
   

    % Boundary Condition for Predicted Velocities
    ut(:,1)=-ut(:,2);
    ut(:,end)=2-ut(:,end-1);
    ut(1,:)=0;
    ut(end,:)=0;
 %   ut(1:iu,ju)=-ut(1:iu,ju+1);
 %   vt(iu,2:ju)=-ut(iu+1,2:ju);
    
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);
    vt(:,1)=0;
    vt(:,end)=0;
 %   vt(iv,1:jv)=-vt(iv+1,1:jv);
 %   vt(2:iv,jv)=-vt(2:iv,jv+1);
    
    % Solving Poisson Equation
    for it=1:maxiter
        for i=2:Nx
            for j=2:Ny
                p(i,j)=w*0.25*(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-(12*dy/(5*dt))*(vt(i,j)-vt(i,j-1))-(12*dx/(5*dt))*(ut(i,j)-ut(i-1,j)))+(1-w)*p(i,j);    
            end
        end
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);

    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    end

    % Velocities after Second RK-3 Step
    u(2:Nx-1,2:Ny)=ut(2:Nx-1,2:Ny)-(dt/(dx*3))*(p(3:Nx,2:Ny)-p(2:Nx-1,2:Ny));

    v(2:Nx,2:Ny-1)=vt(2:Nx,2:Ny-1)-(dt/(dy*3))*(p(2:Nx,3:Ny)-p(2:Nx,2:Ny-1));
    
    %RK-3 Second Step
    
    % Applying Boundary Conditions for Velocites
    u(:,1)=-u(:,2);
    u(:,end)=2-u(:,end-1);
    u(1,:)=0;
    u(end,:)=0;
%    u(1:iu,ju)=-u(1:iu,ju+1);
%    u(iu,2:ju)=-u(iu+1,2:ju);
    
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);
    v(:,1)=0;
    v(:,end)=0;
%    v(iv,1:jv)=-v(iv+1,1:jv);
%    v(2:iv,jv)=-v(2:iv,jv+1);

    % Calculating Predicted Velocities

    ue1=(u(2:Nx-1,2:Ny)+u(3:Nx,2:Ny))/2;
    uw1=(u(2:Nx-1,2:Ny)+u(1:Nx-2,2:Ny))/2;
    us1=(u(2:Nx-1,2:Ny)+u(2:Nx-1,1:Ny-1))/2;
    un1=(u(2:Nx-1,2:Ny)+u(2:Nx-1,3:Ny+1))/2;
    vs1=(v(2:Nx-1,1:Ny-1)+v(3:Nx,1:Ny-1))/2;
    vn1=(v(2:Nx-1,2:Ny)+v(3:Nx,2:Ny))/2;
        
    G1(2:Nx-1,2:Ny)=-(5/9)*G1(2:Nx-1,2:Ny)+(F1C(ue1,uw1,us1,un1,vs1,vn1,dx,dy)+FV(u(2:Nx-1,2:Ny),u(3:Nx,2:Ny),u(1:Nx-2,2:Ny)...
                             ,u(2:Nx-1,3:Ny+1),u(2:Nx-1,1:Ny-1),dx,dy,mu));
    ut(2:Nx-1,2:Ny)=u(2:Nx-1,2:Ny)+(15*dt/16)*G1(2:Nx-1,2:Ny); 
            
    vn12=(v(2:Nx,2:Ny-1)+v(2:Nx,3:Ny))/2;
    vs12=(v(2:Nx,2:Ny-1)+v(2:Nx,1:Ny-2))/2;
    ve12=(v(2:Nx,2:Ny-1)+v(3:Nx+1,2:Ny-1))/2;
    vw12=(v(2:Nx,2:Ny-1)+v(1:Nx-1,2:Ny-1))/2;
    ue12=(u(2:Nx,2:Ny-1)+u(2:Nx,3:Ny))/2;
    uw12=(u(1:Nx-1,2:Ny-1)+u(1:Nx-1,3:Ny))/2;
        
    G2(2:Nx,2:Ny-1)=-(5/9)*G2(2:Nx,2:Ny-1)+(F2C(vn12,vs12,ve12,vw12,ue12,uw12,dx,dy)+FV(v(2:Nx,2:Ny-1),v(3:Nx+1,2:Ny-1),...
                             v(1:Nx-1,2:Ny-1),v(2:Nx,3:Ny),v(2:Nx,1:Ny-2),dx,dy,mu));
    vt(2:Nx,2:Ny-1)=v(2:Nx,2:Ny-1)+(15*dt/16)*G2(2:Nx,2:Ny-1);
    
    % Boundary Condition for Predicted Velocities
    ut(:,1)=-ut(:,2);
    ut(:,end)=2-ut(:,end-1);
    ut(1,:)=0;
    ut(end,:)=0;
 %   ut(1:iu,ju)=-ut(1:iu,ju+1);
 %   vt(iu,2:ju)=-ut(iu+1,2:ju);
    
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);
    vt(:,1)=0;
    vt(:,end)=0;
 %   vt(iv,1:jv)=-vt(iv+1,1:jv);
 %   vt(2:iv,jv)=-vt(2:iv,jv+1);
    
    % Solving Poisson Equation
    for it=1:maxiter
        for i=2:Nx
            for j=2:Ny
                 p(i,j)=w*0.25*(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-(12*dy/(5*dt))*(vt(i,j)-vt(i,j-1))-(12*dx/(5*dt))*(ut(i,j)-ut(i-1,j)))+(1-w)*p(i,j);   
            end
        end
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);

    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    end

    % Velocities after Second RK-3 Step
    u(2:Nx-1,2:Ny)=ut(2:Nx-1,2:Ny)-(5*dt/(dx*12))*(p(3:Nx,2:Ny)-p(2:Nx-1,2:Ny));

    v(2:Nx,2:Ny-1)=vt(2:Nx,2:Ny-1)-(5*dt/(dy*12))*(p(2:Nx,3:Ny)-p(2:Nx,2:Ny-1));

    % RK-3 Third Step
    
    % Applying Boundary Conditions for Velocites
    u(:,1)=-u(:,2);
    u(:,end)=2-u(:,end-1);
    u(1,:)=0;
    u(end,:)=0;
 %   u(1:iu,ju)=-u(1:iu,ju+1);
 %   u(iu,2:ju)=-u(iu+1,2:ju);
    
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);
    v(:,1)=0;
    v(:,end)=0;
 %   v(iv,1:jv)=-v(iv+1,1:jv);
 %   v(2:iv,jv)=-v(2:iv,jv+1);

    % Calculating Predicted Velocities

    ue1=(u(2:Nx-1,2:Ny)+u(3:Nx,2:Ny))/2;
    uw1=(u(2:Nx-1,2:Ny)+u(1:Nx-2,2:Ny))/2;
    us1=(u(2:Nx-1,2:Ny)+u(2:Nx-1,1:Ny-1))/2;
    un1=(u(2:Nx-1,2:Ny)+u(2:Nx-1,3:Ny+1))/2;
    vs1=(v(2:Nx-1,1:Ny-1)+v(3:Nx,1:Ny-1))/2;
    vn1=(v(2:Nx-1,2:Ny)+v(3:Nx,2:Ny))/2;
        
    G1(2:Nx-1,2:Ny)=-(153/128)*G1(2:Nx-1,2:Ny)+(F1C(ue1,uw1,us1,un1,vs1,vn1,dx,dy)+FV(u(2:Nx-1,2:Ny),u(3:Nx,2:Ny),u(1:Nx-2,2:Ny)...
                             ,u(2:Nx-1,3:Ny+1),u(2:Nx-1,1:Ny-1),dx,dy,mu));
    ut(2:Nx-1,2:Ny)=u(2:Nx-1,2:Ny)+(8*dt/15)*G1(2:Nx-1,2:Ny);
   
            
    vn12=(v(2:Nx,2:Ny-1)+v(2:Nx,3:Ny))/2;
    vs12=(v(2:Nx,2:Ny-1)+v(2:Nx,1:Ny-2))/2;
    ve12=(v(2:Nx,2:Ny-1)+v(3:Nx+1,2:Ny-1))/2;
    vw12=(v(2:Nx,2:Ny-1)+v(1:Nx-1,2:Ny-1))/2;
    ue12=(u(2:Nx,2:Ny-1)+u(2:Nx,3:Ny))/2;
    uw12=(u(1:Nx-1,2:Ny-1)+u(1:Nx-1,3:Ny))/2;
        
    G2(2:Nx,2:Ny-1)=-(153/128)*G2(2:Nx,2:Ny-1)+(F2C(vn12,vs12,ve12,vw12,ue12,uw12,dx,dy)+FV(v(2:Nx,2:Ny-1),v(3:Nx+1,2:Ny-1),...
                             v(1:Nx-1,2:Ny-1),v(2:Nx,3:Ny),v(2:Nx,1:Ny-2),dx,dy,mu));
    vt(2:Nx,2:Ny-1)=v(2:Nx,2:Ny-1)+(8*dt/15)*G2(2:Nx,2:Ny-1);

    % Boundary Condition for Predicted Velocities
    ut(:,1)=-ut(:,2);
    ut(:,end)=2-ut(:,end-1);
    ut(1,:)=0;
    ut(end,:)=0;
 %   ut(1:iu,ju)=-ut(1:iu,ju+1);
 %   vt(iu,2:ju)=-ut(iu+1,2:ju);
    
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);
    vt(:,1)=0;
    vt(:,end)=0;
 %   vt(iv,1:jv)=-vt(iv+1,1:jv);
 %   vt(2:iv,jv)=-vt(2:iv,jv+1);

    % Solving Poisson Equation
    for it=1:maxiter
        for i=2:Nx
            for j=2:Ny
                 p(i,j)=w*0.25*(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-(12*dy/(5*dt))*(vt(i,j)-vt(i,j-1))-(12*dx/(5*dt))*(ut(i,j)-ut(i-1,j)))+(1-w)*p(i,j);    
            end
        end
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);

    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    end

    % Velocities after Second RK-3 Step
    u(2:Nx-1,2:Ny)=ut(2:Nx-1,2:Ny)-(dt/(dx*4))*(p(3:Nx,2:Ny)-p(2:Nx-1,2:Ny));

    v(2:Nx,2:Ny-1)=vt(2:Nx,2:Ny-1)-(dt/(dy*4))*(p(2:Nx,3:Ny)-p(2:Nx,2:Ny-1));
    %RK-3 Complete


% Calculating the Velocities and Vorticity at Grid points

uu(1:Nx,1:Ny)=0.5*(u(1:Nx,2:Ny+1)+u(1:Nx,1:Ny));
vv(1:Nx,1:Ny)=0.5*(v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny));
ww(1:Nx,1:Ny)=(u(1:Nx,2:Ny+1)-u(1:Nx,1:Nx)-v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny))/(2*dx);

% Plotting Velocity Vectors and Vorticity Contours
hold off,quiver(flipud(rot90(uu)),flipud(rot90(vv)),'k');hold on;
contour(flipud(rot90(ww)),100),axis equal; axis([1 Nx 1 Ny]);
%fill([xb,xbo],[yb,ybo],'w');
%contour(flipud(rot90(p)),100),axis equal; axis([1 Nx 1 Ny]);
% hold off
% contourf(x,y,flipud(rot90(ww)))
% hold on;
drawnow
pause(0.01)

t=t+dt % Advancing Time Step
end

%% Post Processing

% Velocity Potential and Stream Function Calculation
[phi,psi]=flowfun(uu',vv');

% Plotting Streamfunction Contours
figure
contourf(X,Y,psi); hold on;  axis equal

% Plotting Streamlines and Velocity Vectors
figure
hold on
plot([0,0,Lx,Lx,0],[0,Ly,Ly,0,0],'k');
quiver(X,Y,uu',vv',2,'b')
%contour(X,Y,psi,'k')
%fill([xb,xbo],[yb,ybo],'w');

axis([0 Lx 0 Ly])
axis equal

