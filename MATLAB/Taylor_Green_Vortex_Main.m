% Code to compute the Navier-Stokes Solution for a Taylor Green Vortex
% Written by Akash Dhruv
clear all
close all
clc

Ly=2*pi;
Lx=2*pi;

mu=1;

U=1;
w=1.1;

Nx=50;
Ny=50;

dx=Lx/(Nx-1);
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
xu=zeros(Nx+1,1);
yu=zeros(Ny+1,1);
u=zeros(Nx+1,Ny+1);
ut=zeros(Nx+1,Ny+1);
G1=zeros(Nx,Ny);

% v-velocity points
xv=zeros(Nx+1,1);
yv=zeros(Ny+1,1);
v=zeros(Nx+1,Ny+1);
vt=zeros(Nx+1,Ny+1);
G2=zeros(Nx,Ny);
%
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

for i=2:Nx
    for j=2:Ny
        xu(i)=x(i);
        yu(j)=(y(j-1)+y(j))/2;
    end
end
xu(1)=xu(2)-dx;
yu(1)=yu(2)-dy;

xu(end)=xu(end-1)+dx;
yu(end)=yu(end-1)+dy;

for i=2:Nx
    for j=2:Ny
        xv(i)=(x(i-1)+x(i))/2;
        yv(j)=y(j);
    end
end

xv(1)=xv(2)-dx;
yv(1)=yv(2)-dy;

xv(end)=xv(end-1)+dx;
yv(end)=yv(end-1)+dy;

for i=2:Nx
    for j=2:Ny
        u(i,j)=-cos(xu(i))*sin(yu(j));
    end
end

for i=2:Nx
    for j=2:Ny
        v(i,j)=sin(xv(i))*cos(yv(j));
    end
end

figure
hold on
plot([0,0,Lx,Lx,0],[0,Ly,Ly,0,0])

for i=1:Nx+1
    for j=1:Ny+1
        plot(xp(i),yp(j),'.b')
    end
end
% for i=1:Nx
%     for j=1:1
%         plot(x(i),y(j),'.r')
%         plot(x(i),y(Ny),'.r')
%     end
% end
for i=1:Nx+1
    for j=1:Ny+1
        plot(xu(i),yu(j),'xk')
    end
end
for i=1:Nx+1
    for j=1:Ny+1
        plot(xv(i),yv(j),'*m')
    end
end
axis equal
%%
t=0;
tmax=1;
nt=300;
dt=(tmax-t)/nt;
[X Y]=meshgrid(x,y);

for tstep=1:nt
t=t+dt
u(1,:)=u(end-1,:);
u(end,:)=u(2,:);
u(:,1)=u(:,end-1);
u(:,end)=u(:,2);

v(1,:)=v(end-1,:);
v(end,:)=v(2,:);
v(:,1)=v(:,end-1);
v(:,end)=v(:,2);
for i=2:Nx
    for j=2:Ny
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


for i=2:Nx
    for j=2:Ny
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
ut(1,:)=ut(end-1,:);
ut(end,:)=ut(2,:);
ut(:,1)=ut(:,end-1);
ut(:,end)=ut(:,2);

vt(1,:)=vt(end-1,:);
vt(end,:)=vt(2,:);
vt(:,1)=vt(:,end-1);
vt(:,end)=vt(:,2);

for i=2:Nx
    for j=2:Ny
    p(i,j)=-0.25*(cos(2*xp(i))+cos(2*yp(j)));
    end
end

p(:,1)=p(:,2);
p(:,end)=p(:,end-1);

p(1,:)=p(2,:);
p(end,:)=p(end-1,:);


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


for i=2:Nx
    for j=2:Ny
        u(i,j)=ut(i,j)-(dt/(3*dx))*(p(i+1,j)-p(i,j));
    end
end


for i=2:Nx
    for j=2:Ny
        v(i,j)=vt(i,j)-(dt/(3*dx))*(p(i,j+1)-p(i,j));
    end
end

u(1,:)=u(end-1,:);
u(end,:)=u(2,:);
u(:,1)=u(:,end-1);
u(:,end)=u(:,2);

v(1,:)=v(end-1,:);
v(end,:)=v(2,:);
v(:,1)=v(:,end-1);
v(:,end)=v(:,2);
for i=2:Nx
    for j=2:Ny
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


for i=2:Nx
    for j=2:Ny
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

ut(1,:)=ut(end-1,:);
ut(end,:)=ut(2,:);
ut(:,1)=ut(:,end-1);
ut(:,end)=ut(:,2);

vt(1,:)=vt(end-1,:);
vt(end,:)=vt(2,:);
vt(:,1)=vt(:,end-1);
vt(:,end)=vt(:,2);


maxiter=500;

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


for i=2:Nx
    for j=2:Ny
        u(i,j)=ut(i,j)-(5*dt/(12*dx))*(p(i+1,j)-p(i,j));
    end
end


for i=2:Nx
    for j=2:Ny
        v(i,j)=vt(i,j)-(5*dt/(12*dx))*(p(i,j+1)-p(i,j));
    end
end

u(1,:)=u(end-1,:);
u(end,:)=u(2,:);
u(:,1)=u(:,end-1);
u(:,end)=u(:,2);

v(1,:)=v(end-1,:);
v(end,:)=v(2,:);
v(:,1)=v(:,end-1);
v(:,end)=v(:,2);
for i=2:Nx
    for j=2:Ny
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


for i=2:Nx
    for j=2:Ny
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
ut(1,:)=ut(end-1,:);
ut(end,:)=ut(2,:);
ut(:,1)=ut(:,end-1);
ut(:,end)=ut(:,2);

vt(1,:)=vt(end-1,:);
vt(end,:)=vt(2,:);
vt(:,1)=vt(:,end-1);
vt(:,end)=vt(:,2);


maxiter=500;

for it=1:maxiter

    for i=2:Nx
        for j=2:Ny
        p(i,j)=w*0.25*(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-(4*dx/(dt))*(vt(i,j)-vt(i,j-1)+ut(i,j)-ut(i-1,j)))+(1-w)*p(i,j);  
        end

    end


p(:,1)=p(:,2);
p(:,end)=p(:,end-1);

p(1,:)=p(2,:);
p(end,:)=p(end-1,:);
end


for i=2:Nx
    for j=2:Ny
        u(i,j)=ut(i,j)-(dt/(4*dx))*(p(i+1,j)-p(i,j));
    end
end


for i=2:Nx
    for j=2:Ny
        v(i,j)=vt(i,j)-(dt/(4*dx))*(p(i,j+1)-p(i,j));
    end
end




uu(1:Nx,1:Ny)=0.5*(u(1:Nx,2:Ny+1)+u(1:Nx,1:Ny));
vv(1:Nx,1:Ny)=0.5*(v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny));
ww(1:Nx,1:Ny)=(u(1:Nx,2:Ny+1)-u(1:Nx,1:Nx)-v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny))/(2*dx);

hold off,quiver(flipud(rot90(uu)),flipud(rot90(vv)),'k');hold on;
contour(flipud(rot90(ww)),100),axis equal; axis([1 Nx 1 Ny]);
%contour(flipud(rot90(p)),100),axis equal; axis([1 Nx 1 Ny]);
% hold off, quiver(X,Y,uu',vv');hold on;
% contour(X,Y,ww')
% axis equal
drawnow
pause(0.01)
end
%%
[xnu,ynu]=meshgrid(xu(2:end-1),yu(2:end-1));
figure
quiver(xnu,ynu,u(2:end-1,2:end-1)',v(2:end-1,2:end-1)');
axis equal

figure
hold on
contour(X,Y,ww')
quiver(X,Y,uu',vv','k')

figure
contour(X,Y,uu')
title('u-numerical')

figure
contour(X,Y,vv')
title('v-numerical')

% [Xv Yv]=meshgrid(xv(2:end-1),yv);
% figure
% contour(Xv,Yv,v(2:end-1,:)')
% 
% [Xp Yp]=meshgrid(xp(2:end-1),yp(2:end-1));
% figure
% contour(Xp,Yp,p(2:end-1,2:end-1)');

for i =1:Nx
    for j=1:Ny
        uann(i,j)=-exp(-2*t)*cos(X(i,j))*sin(Y(i,j));
    end
end

for i =1:Nx
    for j=1:Ny
        vann(i,j)=exp(-2*t)*sin(X(i,j))*cos(Y(i,j));
    end
end

[Xp,Yp]=meshgrid(xp,yp);
for i =1:Nx+1
    for j=1:Ny+1
        pann(i,j)=-0.25*exp(-4*t)*(cos(2*Xp(i,j))+cos(2*Yp(i,j)));
    end
end



figure
contour(X,Y,uann)
title('u-analytical')

figure
contour(X,Y,vann)
title('v-analytical')

figure
contour(Xp,Yp,pann)
title('p-analytical')

figure
contour(Xp,Yp,p');
title('p-numerical')

diffu=uu'-uann;

diffv=vv'-vann;

diffp=p'-pann;

figure
mesh(X,Y,diffu)
title('U_error')
axis equal

figure
mesh(X,Y,diffv)
title('V_error')
axis equal

figure
mesh(Xp,Yp,diffp)
title('P_error')
axis equal
%[startx,starty]=meshgrid(startx,starty);
% [phi,psi]=flowfun(uu',vv');
% [Xu Yu]=meshgrid(xu,yu(2:end-1));
% figure
% contour(Xu,Yu,u(:,2:end-1)'); hold on; axis equal
% 
% % figure
% % hold on
% % plot([0,0,1.5,1.5,0],[0,1.5,1.5,0,0],'k');
% % quiver(X,Y,uu',vv',2,'b')
% % % streamline(X,Y,uu',vv',startx,starty)
% % contour(X,Y,psi,'k')
% % axis([0 2*pi 0 2*pi])
% % axis equal
% 
% % figure
% % mesh(xp(2:Nx),yp(2:Nx),p(2:Nx,2:Ny))
