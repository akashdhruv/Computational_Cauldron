% Code to compute the Navier-Stokes Solution for a Lid-Driven Cavity
% Written by Akash Dhruv
clear all
close all
clc

Ly=1; % Domain Lengths
Lx=1;

r=0.05;

%mu=20.94*(10^-6); %m^2/s

mu=0.001;

rho=1; % kg/m^3
Re=300;
L=2*r;

% U=(Re*mu)/L;
U=1;

w=1.1; % SOR Relaxation Factor

Nx=150; % Grid Size
Ny=150;

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
% Plotting the Grid


dis=linspace(0,pi,20);
xcircle=0.3+r*cos(dis);
ycircle=sqrt(r^2-(xcircle-0.3).^2);

xcircle=[xcircle,fliplr(xcircle(1:end-1))];
ycircle=[ycircle,-fliplr(ycircle(1:end-1))];

ycircle=ycircle+0.5;

plot(xcircle,ycircle,'-k.')

egu=1;
for i=1:length(xu)
    for j=1:length(yu)
        
        if(sqrt((xu(i)-0.3)^2+(yu(j)-0.5)^2)<=r)
            iu(egu,:)=[i,j];
            egu=egu+1;
        end
        
    end
end
% for i=1:egu-1
%     plot(xu(iu(i,1)),yu(iu(i,2)),'.b')
% end
egv=1;
for i=1:length(xv)
    for j=1:length(yv)
        
        if(sqrt((xv(i)-0.3)^2+(yv(j)-0.5)^2)<=r)
            iv(egv,:)=[i,j];
            egv=egv+1;
        end
        
    end
end
% for i=1:egv-1
%     plot(xv(iv(i,1)),yv(iv(i,2)),'xk')
% end
% %%
% bul=find(iu(:,1)==min(iu(:,1)));
% bur=find(iu(:,1)==max(iu(:,1)));
% bub=find(iu(:,2)==min(iu(:,2)));
% but=find(iu(:,2)==max(iu(:,2)));
% 
% bvl=find(iv(:,1)==min(iv(:,1)));
% bvr=find(iv(:,1)==max(iv(:,1)));
% bvb=find(iv(:,2)==min(iv(:,2)));
% bvt=find(iv(:,2)==max(iv(:,2)));
% 
% bu=[iu(bul,1)-1,iu(bul,2);iu(bur,1)+1,iu(bur,2);iu(bub,1),iu(bub,2)-1;iu(but,1),iu(but,2)+1];
% 
% for i=1:length(bu)
%     plot(xu(bu(i,1)),yu(bu(i,2)),'.b')
% end
% 
% bv=[iv(bvl,1)-1,iv(bvl,2);iv(bvr,1)+1,iv(bvr,2);iv(bvb,1),iv(bvb,2)-1;iv(bvt,1),iv(bvt,2)+1];
% 
% % for i=1:length(bv)
% %     plot(xv(bv(i,1)),yv(bv(i,2)),'xk')
% % end
% %%

bul=min(iu(:,1));
bur=max(iu(:,1));
bub=min(iu(:,2));
but=max(iu(:,2));

cux=linspace(bul,bur,bur-bul+1);

count=1;
for i=1:length(cux)
bu(count,:)=[cux(i),min(iu(find(iu(:,1)==cux(i)),2))-1];
count=count+1;
bu(count,:)=[cux(i),max(iu(find(iu(:,1)==cux(i)),2))+1];
count=count+1;
end
 

cuy=linspace(bub,but,but-bub+1);

for i=1:length(cuy)
bu(count,:)=[min(iu(find(iu(:,2)==cuy(i)),1))-1,cuy(i)];
count=count+1;
bu(count,:)=[max(iu(find(iu(:,2)==cuy(i)),1))+1,cuy(i)];
count=count+1;
end
    
for i=1:length(bu)
    plot(xu(bu(i,1)),yu(bu(i,2)),'.b')
end    
    
bvl=min(iv(:,1));
bvr=max(iv(:,1));
bvb=min(iv(:,2));
bvt=max(iv(:,2));

cvx=linspace(bvl,bvr,bvr-bvl+1);

covnt=1;
for i=1:length(cvx)
bv(covnt,:)=[cvx(i),min(iv(find(iv(:,1)==cvx(i)),2))-1];
covnt=covnt+1;
bv(covnt,:)=[cvx(i),max(iv(find(iv(:,1)==cvx(i)),2))+1];
covnt=covnt+1;
end
 

cvy=linspace(bvb,bvt,bvt-bvb+1);

for i=1:length(cvy)
bv(covnt,:)=[min(iv(find(iv(:,2)==cvy(i)),1))-1,cvy(i)];
covnt=covnt+1;
bv(covnt,:)=[max(iv(find(iv(:,2)==cvy(i)),1))+1,cvy(i)];
covnt=covnt+1;
end
    
for i=1:length(bv)
    plot(xv(bv(i,1)),yv(bv(i,2)),'xk')
end    


axis equal
%%

t=0; % t-initial
tmax=5; % t-final
dt=0.45*min(dx,dy);
nt=ceil((tmax-t)/dt);
maxiter=100;
[X Y]=meshgrid(x,y);

figure
% Starting time loop
for tstep=1:nt
     
    % RK3 First Step
    
     % Applying Boundary Conditions for Velocites

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
    
    
    for i=1:length(bu)
        ut(bu(i,1),bu(i,2))=0;
    end
    for i=1:length(bv)
        vt(bv(i,1),bv(i,2))=0;
    end

    % Boundary Condition for Predicted Velocities
    ut(:,1)=2*U-ut(:,end-1);
    ut(:,end)=2*U-ut(:,2);
    ut(1,2:end-1)=U;
    ut(end,2:end-1)=ut(end,2:end-1)-U*(ut(end,2:end-1)-ut(end-1,2:end-1))*(dt/dx);
    
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);
    vt(2:end-1,1)=0;
    vt(2:end-1,end)=0;

%     ut(:,1)=-ut(:,2);
%     ut(:,end)=2*U-ut(:,end-1);
%     
%     
%     vt(1,:)=0-vt(2,:);
%     vt(end,:)=0-vt(end-1,:);
%     vt(2:end-1,end)=- 0.5;
    
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
    u(:,1)=2*U-u(:,end-1);
    u(:,end)=2*U-u(:,2);
    u(1,2:end-1)=U;
    u(end,2:end-1)=u(end,2:end-1)-U*(u(end,2:end-1)-u(end-1,2:end-1))*(dt/dx);
    
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);
    v(2:end-1,1)=0;
    v(2:end-1,end)=0;
    
%     u(:,1)=-u(:,2);
%     u(:,end)=2*U-u(:,end-1);
%     
%     
%     v(1,:)=0-v(2,:);
%     v(end,:)=0-v(end-1,:);
%     v(2:end-1,end)= -0.5;



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
    
    
    for i=1:length(bu)
        ut(bu(i,1),bu(i,2))=0;
    end
    for i=1:length(bv)
        vt(bv(i,1),bv(i,2))=0;
    end
    
    % Boundary Condition for Predicted Velocities
    ut(:,1)=2*U-ut(:,end-1);
    ut(:,end)=2*U-ut(:,2);
    ut(1,2:end-1)=U;
    ut(end,2:end-1)=ut(end,2:end-1)-U*(ut(end,2:end-1)-ut(end-1,2:end-1))*(dt/dx);
    
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);
    vt(2:end-1,1)=0;
    vt(2:end-1,end)=0;

%     ut(:,1)=-ut(:,2);
%     ut(:,end)=2*U-ut(:,end-1);
%     
%     
%     vt(1,:)=0-vt(2,:);
%     vt(end,:)=0-vt(end-1,:);
%     vt(2:end-1,end)= -0.5;


    
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
    u(:,1)=2*U-u(:,end-1);
    u(:,end)=2*U-u(:,2);
    u(1,2:end-1)=U;
    u(end,2:end-1)=u(end,2:end-1)-U*(u(end,2:end-1)-u(end-1,2:end-1))*(dt/dx);
    
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);
    v(2:end-1,1)=0;
    v(2:end-1,end)=0;

%     u(:,1)=-u(:,2);
%     u(:,end)=2*U-u(:,end-1);
%     
%     
%     v(1,:)=0-v(2,:);
%     v(end,:)=0-v(end-1,:);
%     v(2:end-1,end)= -0.5;


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
     
    
    for i=1:length(bu)
        ut(bu(i,1),bu(i,2))=0;
    end
    for i=1:length(bv)
        vt(bv(i,1),bv(i,2))=0;
    end
    
    % Boundary Condition for Predicted Velocities
    ut(:,1)=2*U-ut(:,end-1);
    ut(:,end)=2*U-ut(:,2);
    ut(1,2:end-1)=U;
    ut(end,2:end-1)=ut(end,2:end-1)-U*(ut(end,2:end-1)-ut(end-1,2:end-1))*(dt/dx);
    
    
    vt(1,:)=0-vt(2,:);
    vt(end,:)=0-vt(end-1,:);
    vt(2:end-1,1)=0;
    vt(2:end-1,end)=0;

%     ut(:,1)=-ut(:,2);
%     ut(:,end)=2*U-ut(:,end-1);
%     
%     
%     vt(1,:)=0-vt(2,:);
%     vt(end,:)=0-vt(end-1,:);
%     vt(2:end-1,end)= -0.5;



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
    
    u(:,1)=2*U-u(:,end-1);
    u(:,end)=2*U-u(:,2);
    u(1,2:end-1)=U;
    u(end,2:end-1)=u(end,2:end-1)-U*(u(end,2:end-1)-u(end-1,2:end-1))*(dt/dx);
    
    
    v(1,:)=0-v(2,:);
    v(end,:)=0-v(end-1,:);
    v(2:end-1,1)=0;
    v(2:end-1,end)=0;


%     u(:,1)=-u(:,2);
%     u(:,end)=2*U-u(:,end-1);
%     
%     
%     v(1,:)=0-v(2,:);
%     v(end,:)=0-v(end-1,:);
%     v(2:end-1,end)= -0.5;

    
    %RK-3 Complete


% Calculating the Velocities and Vorticity at Grid points

uu(1:Nx,1:Ny)=0.5*(u(1:Nx,2:Ny+1)+u(1:Nx,1:Ny));
vv(1:Nx,1:Ny)=0.5*(v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny));
ww(1:Nx,1:Ny)=(u(1:Nx,2:Ny+1)-u(1:Nx,1:Nx)-v(2:Nx+1,1:Ny)+v(1:Nx,1:Ny))/(2*dx);



t=t+dt % Advancing Time Step

% Plotting Velocity Vectors and Vorticity Contours
hold off
hold on;
[h h]=contourf(X,Y,ww',200);
set(h,'LineStyle','none');
fill(xcircle,ycircle,'w')
axis equal
%fill([xb,xbo],[yb,ybo],'w');
%contour(flipud(rot90(p)),100),axis equal; axis([1 Nx 1 Ny]);
% hold off
% contourf(x,y,flipud(rot90(ww)))
% hold on;
drawnow
pause(0.01)
   
   
   
   
end

%% Post Processing
figure
hold on
axis([0 Lx 0 Ly])
[h h]=contourf(X,Y,sqrt(uu.^2+vv.^2)',200);
set(h,'LineStyle','none');
quiver(X,Y,uu',vv',2,'k')
fill(xcircle,ycircle,'w');
axis equal

figure
hold on
axis([0 Lx 0 Ly])
[h h]=contourf(X,Y,ww',200);
set(h,'LineStyle','none');
% quiver(X,Y,uu',vv',2,'k')
fill(xcircle,ycircle,'w');
axis equal



