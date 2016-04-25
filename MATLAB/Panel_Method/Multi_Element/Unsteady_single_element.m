clear all
close all
clc

% Defining Geometry and Freestream Conditions
M=20;
alpha=10*pi/180;
Q_inf=1;
U_inf=Q_inf*cos(alpha);
W_inf=Q_inf*sin(alpha);

[xu,zu,xl,zl]= naca4m('4412',1,M);
xm=[xl; flipud(xu(1:end-1))];
zm=[zl; flipud(zu(1:end-1))];
zm(end)=0;
figure;
plot(xm,zm,'b.-');
axis equal;

N=length(xm);

% Setting the Number of time steps Nt and time t
Nt=200;
t=10;
dt=t/Nt;

% Getting collocation points
xc=(xm(1:end-1)+xm(2:end))/2;
zc=(zm(1:end-1)+zm(2:end))/2;

% Setting vectors for wake geometry and and wake doublet strengths
% Intial zeros because no wake is shed
MuW=zeros(Nt-1,1);
xw=zeros(Nt,1);
zw=zeros(Nt,1);

% Number of Collcation points = number of panels-1
n_main=N-1;

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
        if i==j;
            C(i,j)= 0.5;
        else
           C(i,j)= PHICD(1,xc(i),zc(i),xm(j),zm(j),xm(j+1),zm(j+1));
        end
    end
end

% End of calculation of Matrix C
% Time loop for wake shedding 'tik' is time counter.
 for tik=1:Nt
 % Getting vector for RHS, for first time step no wake is shed so
 % contribution of wakes is zero
 RHS=zeros(n_main,1);
  for i=1:N-1;
      q=0;
      % Loop to get wake influences
      for j=1:Nt-1
          q=q+PHICD(MuW(j),xc(i),zc(i),xw(j),zw(j),xw(j+1),zw(j+1));
      end
      RHS(i,1)=-[xc(i) zc(i)]*[U_inf;W_inf]-q;
  end
  
% Getting Doublet Strengths for body panels  
Mu=C\RHS;


% Calculating Cp
Cp=zeros(N-1,1); 
for i = 1:N-1

	if i == 1
		R = sqrt((xc(2)-xc(1))^2+(zc(2)-zc(1))^2);
		v_loc_t = (Mu(2)-Mu(1))/R;
	elseif i == N-1
		R = sqrt((xc(N-1)-xc(N-2))^2+(zc(N-1)-zc(N-2))^2);
		v_loc_t = (Mu(N-1)-Mu(N-2))/R;
	else
		R = sqrt((xc(i+1)-xc(i-1))^2+(zc(i+1)-zc(i-1))^2);
		v_loc_t = (Mu(i+1)-Mu(i-1))/R;
	end
	v_t =v_loc_t;
	Cp(i) = 1 - ((v_t)^2);
end



% title('Cp');
% plot(xc(1:end),Cp,'-r');
% set(gca,'YDir','reverse');
% ylabel('C_p');
% xlabel('x');
% drawnow

% Getting the Wake strengths at each loop by enforcing Kutta Condition
% mu_wake=mu_ut-mu_lt

if(tik==1)
    MuW(1)=Mu(end)-Mu(1);
else
    MuW(2:end)=MuW(1:end-1);
    MuW(1)=Mu(end)-Mu(1);
end

% Shedding Wake panels at each time step
if(tik==1)
    xw(tik)=xm(end);
    zw(tik)=zm(end);
    xw(tik+1)=xm(end)+tik*Q_inf*dt;
    zw(tik+1)=zm(end);
else
    xw(1:tik)=xw(1:tik)+Q_inf*dt;
    zw(1:tik)=zw(1:tik);
    xw(2:end)=xw(1:end-1);
    zw(2:end)=zw(1:end-1);
    xw(1)=xm(end);
    zw(1)=zm(end);
end
if(tik>1)
 
 % Calculating Induced Velocities on Wake panels and Wake shapes at each
 % time step
 uw=zeros(tik-1,1);
 ww=zeros(tik-1,1);
 xcw=(xw(1:tik-1)+xw(2:tik))/2;
 zcw=(zw(1:tik-1)+zw(2:tik))/2;
 for i=1:tik-1;
    u=0; w=0;
    for j=1:max(n_main,Nt);
        if(j<=n_main)
        [up,wp]=doublet_2Dc(Mu(j),xcw(i),zcw(i),xm(j),zm(j),xm(j+1),zm(j+1));
        else
         up=0; wp=0;
        end
        if (j<=Nt-1)
         [upw,wpw]=doublet_2Dc(MuW(j),xcw(i),zcw(i),xw(j),zw(j),xw(j+1),zw(j+1));
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
for i=2:tik
    zw(i)=zw(i-1)+(xw(i)-xw(i-1))*(ww(i-1)/uw(i-1));
    xw(i)=xw(i-1)+(zw(i)-zw(i-1))*(uw(i-1)/ww(i-1));
end

% Plotting the results
plot(xm,zm,'b',xw(1:tik),zw(1:tik),'k');
axis equal;
drawnow
end
end

