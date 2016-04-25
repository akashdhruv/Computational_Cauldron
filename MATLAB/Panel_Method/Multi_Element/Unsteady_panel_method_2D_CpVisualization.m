clear all
close all
clc


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


Nt=100;
t=20;
dt=t/Nt;

xc=(xm(1:end-1)+xm(2:end))/2;
zc=(zm(1:end-1)+zm(2:end))/2;

MuW=zeros(Nt-1,1);
xw=zeros(Nt,1);
zw=zeros(Nt,1);

n_main=N-1;

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


 for tik=1:Nt
 RHS=zeros(n_main,1);
  for i=1:N-1;
      q=0;
      for j=1:Nt-1
          q=q+PHICD(MuW(j),xc(i),zc(i),xw(j),zw(j),xw(j+1),zw(j+1));
      end
      RHS(i,1)=-[xc(i) zc(i)]*[U_inf;W_inf]-q;
  end
  
Mu=C\RHS;

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



title('Cp');
plot(xc,Cp,'-r');
set(gca,'YDir','reverse');
ylabel('C_p');
xlabel('x');
drawnow

if(tik==1)
    MuW(1)=Mu(end)-Mu(1);
else
    MuW(2:end)=MuW(1:end-1);
    MuW(1)=Mu(end)-Mu(1);
end

if(tik==1)
    xw(tik)=xm(end);
    zw(tik)=zm(end);
end
xw(tik+1)=xm(end)+tik*Q_inf*dt;
zw(tik+1)=zm(end);
 end
 uw=zeros(Nt-1,1);
 ww=zeros(Nt-1,1);
xcw=(xw(1:end-1)+xw(2:end))/2;
zcw=(zw(1:end-1)+zw(2:end))/2;
MuW=[Mu(end)-Mu(1);MuW];
 for i=1:Nt;
    u=0; w=0;
    for j=1:max(n_main,Nt);
        if(j<=n_main)
        [up,wp]=doublet_2Dc(Mu(j),xcw(i),zcw(i),xm(j),zm(j),xm(j+1),zm(j+1));
        else
         up=0; wp=0;
        end
        if (j<=Nt)
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
for i=2:Nt+1
    zw(i)=zw(i-1)+(xw(i)-xw(i-1))*(ww(i-1)/uw(i-1));
    xw(i)=xw(i-1)+(zw(i)-zw(i-1))*(uw(i-1)/ww(i-1));
end
figure;
plot(xm,zm,'b',xw,zw,'k');
axis equal;

