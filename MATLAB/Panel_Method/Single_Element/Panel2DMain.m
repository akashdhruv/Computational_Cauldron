% Sample program to calculate Coefficient of Pressure for a two dimensional
% airfoil using a constant strength doublet panel method with dirichlet
% boundary condition

% Written by Akash Dhruv

% Distributed under Creative Commons License 

clear all
close all
clc

M=30;
N=2*M;

%%%%%%%%%%%%%%%%%%%%% DEFINE AIRFOL and ANGLE OF ATTACK %%%%%%%%%%%%%%%%
[xu,zu,xl,zl]= naca4m('4412',1,M);
alpha = 5; % Degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:M+1;
for i=1:M;
    if  xu(i+1)<xu(i); %sorting the array for x cordinates of upper surface
        xu(i)=xu(i)+xu(i+1);% will incorporate the flipup and fliplr next time
        xu(i+1)= xu(i)-xu(i+1);
        xu(i)=xu(i)-xu(i+1);
        zu(i)=zu(i)+zu(i+1);
        zu(i+1)= zu(i)-zu(i+1);
        zu(i)=zu(i)-zu(i+1);
    end
end
end
xplot= [xl(1:M);xu]; %arranging panels from lowerside of trailing edge to upper side of trailing edge
zplot= [zl(1:M);zu];
zplot(end)=0;
figure;
plot(xplot,zplot,'b.-'); % plotting the airfoil
axis equal;

%Calculating collocation points
    xmid=(xplot(1:end-1)+xplot(2:end))/2;
    zmid=(zplot(1:end-1)+zplot(2:end))/2;
    
  n=zeros(2,N);
  C= zeros(N+1,N+1); % Matrix for influence coefficeints
  for i=1:N; %calculating influence coefficients Cij
      al(i)= atan2(zplot(i+1)-zplot(i),xplot(i+1)-xplot(i));
      n(:,i)= [sin(al(i)) cos(al(i))];
   
      for j=1:N;
          if i==j;
            C(i,j)= 0.5;% influence coefficent when i=j
          else
           C(i,j)= PHICD(1,xmid(i),zmid(i),xplot(j),zplot(j),xplot(j+1),zplot(j+1)); %calling function PHICD 
          end
      end
      %Influence due to wake panels
      C(i,N+1)= PHICD(1,xmid(i),zmid(i),xplot(1),zplot(1),10966,zplot(1));
  end
   
  %Explicit Kutta Condition
  C(N+1,1)=1;
  C(N+1,N)=-1;
  C(N+1,N+1)=1;


  % Free stream velocity comoponents
  alpha=alpha*pi/180;
  Q_inf=1;
  U_inf=Q_inf*cos(alpha);
  V_inf=Q_inf*sin(alpha);

  RHS=zeros(N+1,1);
  for i=1:N; % calculating RHS, i.e. potential at collocation points
      RHS(i,1)= -[xmid(i) zmid(i)]*[U_inf;V_inf];
  end
  
  DoubS=C\RHS; % Doublet Strength

 
for i = 1:N
% 	[u,w] = doublet_2Dc(DoubS(N+1),xmid(i),zmid(i),xplot(1),zplot(1),inf,inf);
% 	v_t = [u w]*[n(2,i);n(1,i)];
% 	for j = 1:N
% 		if j ~= i
% 			[u,w] = doublet_2Dc(DoubS(j),xmid(i),zmid(i),xplot(j),zplot(j),xplot(j+1),zplot(j+1));
% 			v_t = v_t + [u w]*[n(2,i);n(1,i)];
% 		end
% 	end
	if i == 1
		R = sqrt((xmid(2)-xmid(1))^2+(zmid(2)-zmid(1))^2);
		v_loc_t = (DoubS(2)-DoubS(1))/R;
	elseif i == N
		R = sqrt((xmid(N)-xmid(N-1))^2+(zmid(N)-zmid(N-1))^2);
		v_loc_t = (DoubS(N)-DoubS(N-1))/R;
	else
		R = sqrt((xmid(i+1)-xmid(i-1))^2+(zmid(i+1)-zmid(i-1))^2);
		v_loc_t = (DoubS(i+1)-DoubS(i-1))/R;
	end
	v_t =v_loc_t;
	%V_inf_t = [U_inf W_inf]*t_vec(:,i);
	Cp(i) = 1 - ((v_t)^2);
end

ind = length(xmid)/2;

figure;
hold on
%title('Angle of Attack al');
plot(xmid(1:ind),Cp(1:ind),'.-k',xmid(ind+1:end),Cp(ind+1:end),'.-r');
set(gca,'YDir','reverse');
ylabel('C_p');
xlabel('x');
legend('lower','upper');
hold off
