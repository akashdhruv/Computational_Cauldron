clear all
close all
clc

% Defining Geometry and Freestream Conditions
M=20;
alpha=10*pi/180;
delta=21*pi/180;
x_disp=1.02;
z_disp=-0.05;

Q_inf=1;
U_inf=Q_inf*cos(alpha);
W_inf=Q_inf*sin(alpha);

[xu,zu,xl,zl]= naca4m('4412',1,M);
xmain=[xl; flipud(xu(1:end-1))];
zmain=[zl; flipud(zu(1:end-1))];
zmain(end)=0;

[xu,zu,xl,zl]  = naca4m('4415',0.5,M);
xflap = [xl; flipud(xu(1:end-1))]; 
zflap = [zl; flipud(zu(1:end-1))];
zflap(end)=0;

x1 = zeros(size(xflap)); 
z1 = zeros(size(zflap)); 

for i = 1:length(xflap) 
    r = sqrt(xflap(i)^2 + zflap(i)^2); 
    theta_r = atan2(zflap(i),xflap(i)); 
    theta = theta_r- delta; 
    x1(i) = r*cos(theta); 
    z1(i) = r*sin(theta);
end

xflap = x1; 
zflap = z1;

xflap=xflap+x_disp;
zflap=zflap+z_disp;

xm=zeros(length(xmain),2);
zm=zeros(length(xmain),2);

xm(:,1)=xmain; zm(:,1)=zmain;
xm(:,2)=xflap; zm(:,2)=zflap;


figure;
plot(xm(:,1),zm(:,1),'b',xm(:,2),zm(:,2),'k');
axis equal;

N=length(xm);

% Setting the Number of time steps Nt and time t
Nt=100;
t=10;
dt=t/Nt;

% Getting collocation points
xc=(xm(1:end-1,:)+xm(2:end,:))/2;
zc=(zm(1:end-1,:)+zm(2:end,:))/2;

figure;
hold on
plot(xm(:,1),zm(:,1),'b',xm(:,2),zm(:,2),'k');
plot(xc(:,1),zc(:,1),'.k',xc(:,2),zc(:,2),'.b');
axis equal;
%
% Setting vectors for wake geometry and and wake doublet strengths
% Intial zeros because no wake is shed
MuW=zeros(Nt-1,2);
xw=zeros(Nt,2);
zw=zeros(Nt,2);

% Number of Collcation points = number of panels-1
n_main=2*(N-1);

% Calculating Influences of Body panels on itself. Initially no wake is
% shed. We need to calculate influence matrix just once. Influence of wakes
% will be added to RHS, kutta condition will be satisfied in calclulating
% wake strengths.

C=zeros(n_main);
al=zeros(n_main);
n=zeros(2,n_main);

for i=1:n_main
%     al(i)= atan2(zm(i+1)-zm(i),xm(i+1)-xm(i));
%     n(:,i)= [sin(al(i)) cos(al(i))];
%     dc=sqrt((zm(i+1)-zm(i)^2)-(xm(i+1)-xm(i))^2);
    for j=1:n_main
        if j==i;
            C(i,j)= 0.5;
        elseif j<=N-1
           C(i,j)= PHICD(1,xc(i),zc(i),xm(j),zm(j),xm(j+1),zm(j+1));
        else
           C(i,j)= PHICD(1,xc(i),zc(i),xm(j+1),zm(j+1),xm(j+2),zm(j+2));    
        end
    end
end

figure;
% End of calculation of Matrix C
% Time loop for wake shedding 'tik' is time counter.
 for tik=1:Nt
 % Getting vector for RHS, for first time step no wake is shed so
 % contribution of wakes is zero
 RHS=zeros(n_main,1);
  for i=1:n_main;
      q=0;
      % Loop to get wake influences
      for j=1:2*(Nt-1)
          if(j<=Nt-1)
          q=q+PHICD(MuW(j),xc(i),zc(i),xw(j),zw(j),xw(j+1),zw(j+1));
          else
          q=q+PHICD(MuW(j),xc(i),zc(i),xw(j+1),zw(j+1),xw(j+2),zw(j+2));
          end
      end
      RHS(i,1)=-[xc(i) zc(i)]*[U_inf;W_inf]-q;
  end
  
% Getting Doublet Strengths for body panels  
Mu=C\RHS;
MU=zeros(N-1,2);
k=0;
for j=1:2
    for i=1:N-1
        k=k+1;
        MU(i,j)=Mu(k);
    end
end


% Calculating Cp
Cp=zeros((N-1),2); 
for i = 1:N-1
    for j=1:2

	if i == 1
		R = sqrt((xc(2,j)-xc(1,j))^2+(zc(2,j)-zc(1,j))^2);
		v_loc_t = (MU(2,j)-MU(1,j))/R;
	elseif i == N-1
		R = sqrt((xc(N-1,j)-xc(N-2,j))^2+(zc(N-1,j)-zc(N-2,j))^2);
		v_loc_t = (MU(N-1,j)-MU(N-2,j))/R;
	else
		R = sqrt((xc(i+1,j)-xc(i-1,j))^2+(zc(i+1,j)-zc(i-1,j))^2);
		v_loc_t = (MU(i+1,j)-MU(i-1,j))/R;
	end
	v_t =v_loc_t;
	Cp(i,j) = 1 - ((v_t)^2);
    end
end

% 
% 
% title('Cp');
% plot(xc(:,1),Cp(:,1),'-r',xc(:,2),Cp(:,2),'-k');
% set(gca,'YDir','reverse');
% ylabel('C_p');
% xlabel('x');
% drawnow

% Getting the Wake strengths at each loop by enforcing Kutta Condition
% mu_wake=mu_ut-mu_lt

if(tik==1)
    MuW(1,1)=MU(end,1)-MU(1,1);
    MuW(1,2)=MU(end,2)-MU(1,2);
else
    MuW(2:end,1)=MuW(1:end-1,1);
    MuW(1,1)=MU(end,1)-MU(1,1);
    
    MuW(2:end,2)=MuW(1:end-1,2);
    MuW(1,2)=MU(end,2)-MU(1,2);
end

% Shedding Wake panels at each time step
if(tik==1)
    xw(tik,1)=xm(end,1);
    zw(tik,1)=zm(end,1);
    xw(tik+1,1)=xm(end,1)+tik*Q_inf*dt;
    zw(tik+1,1)=zm(end,1);
    
    xw(tik,2)=xm(end,2);
    zw(tik,2)=zm(end,2);
    xw(tik+1,2)=xm(end,2)+tik*Q_inf*dt;
    zw(tik+1,2)=zm(end,2);
else
    xw(1:tik,1)=xw(1:tik,1)+Q_inf*dt;
    zw(1:tik,1)=zw(1:tik,1)+Q_inf*dt*tan(alpha);
    xw(2:end,1)=xw(1:end-1,1);
    zw(2:end,1)=zw(1:end-1,1);
    xw(1,1)=xm(end,1);
    zw(1,1)=zm(end,1);
    
    xw(1:tik,2)=xw(1:tik,2)+Q_inf*dt;
    zw(1:tik,2)=zw(1:tik,2)+Q_inf*dt*tan(alpha);
    xw(2:end,2)=xw(1:end-1,2);
    zw(2:end,2)=zw(1:end-1,2);
    xw(1,2)=xm(end,2);
    zw(1,2)=zm(end,2);
end
if(tik>1)
 
 % Calculating Induced Velocities on Wake panels and Wake shapes at each
 % time step
 xcw=(xw(1:tik-1,:)+xw(2:tik,:))/2;
 zcw=(zw(1:tik-1,:)+zw(2:tik,:))/2;
 uw=zeros(2*(tik-1),1);
 ww=zeros(2*(tik-1),1);
 for i=1:2*(tik-1);
    u=0; w=0;
    for j=1:max(n_main,2*(Nt-1));
        
        if(j<=N-1)
        [up,wp]=doublet_2Dc(Mu(j),xcw(i),zcw(i),xm(j),zm(j),xm(j+1),zm(j+1));
        elseif(j>N-1 && j<=n_main)
        [up,wp]=doublet_2Dc(Mu(j),xcw(i),zcw(i),xm(j+1),zm(j+1),xm(j+2),zm(j+2));
        else
         up=0; wp=0;
        end
        
        if (j<=Nt-1)
         [upw,wpw]=doublet_2Dc(MuW(j),xcw(i),zcw(i),xw(j),zw(j),xw(j+1),zw(j+1));
        elseif(j>Nt-1 && j<=2*(Nt-1))
         [upw,wpw]=doublet_2Dc(MuW(j),xcw(i),zcw(i),xw(j+1),zw(j+1),xw(j+2),zw(j+2)); 
        else
         upw=0;
         wpw=0;
        end
         
      u=u-up-upw;
      w=w-wp-wpw;
    end
    uw(i)=u+U_inf;
    ww(i)=w+W_inf;
 end
 

k=0;
for j=1:2
    for i=2:tik
        k=k+1;
    zw(i,j)=zw(i-1,j)+(xw(i,j)-xw(i-1,j))*(ww(k)/uw(k));
    xw(i,j)=xw(i-1,j)+(zw(i,j)-zw(i-1,j))*(uw(k)/ww(k));
    end
end
% 
% Plotting the results
plot(xm(:,1),zm(:,1),'b',xw(1:tik,1),zw(1:tik,1),'k',xm(:,2),zm(:,2),'b',xw(1:tik,2),zw(1:tik,2),'k');
axis equal;
drawnow
end
end
figure
title('Cp');
plot(xc(:,1),Cp(:,1),'-r',xc(:,2),Cp(:,2),'-k');
set(gca,'YDir','reverse');
ylabel('C_p');
xlabel('x');
