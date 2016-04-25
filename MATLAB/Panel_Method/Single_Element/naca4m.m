function [xu,zu,xl,zl]  = naca4m(NACA,c,n_panels)
%
% [xu,zu,xl,zl]  = naca4(NACA,c,n_panels)
%
% [Inputs]
% NACA = 4-digit NACA airfoil number (input as a string)
% c = chord
% n_panels = number of panels on each surface
%
% [Outputs]
% xu = x-coordinates of upper surface
% zu = z-coordinates of upper surface
% xl = x-coordinates of lower surface
% zl = z-coordinates of lower surface
%
% Written by Adam Wickenheiser
%



% decode NACA airfoil number
H = str2double(NACA(1))/100;
p = str2double(NACA(2))/10;
T = str2double(NACA(3:4))/100;

% x distribution
beta = (0:(pi/n_panels):pi)';
xc = c*(1-.5*(1-cos(beta)));

% thickness distribution
thdis = 5*T*c*(0.2969*sqrt(xc/c)-0.126*xc/c-0.3537*(xc/c).^2 +0.2843*(xc/c).^3-0.1015*(xc/c).^4);

% camberline
if p ~= 0 && H ~= 0 
 I1 = find(xc <= p*c);
 I2 = find(xc > p*c);
 camberline(I1,1) = (H/p^2)*xc(I1).*(2*p-xc(I1)/c);
 camberline(I2,1) = (H/(1-p)^2)*(c-xc(I2)).*(1+xc(I2)/c-2*p);
end

% airfoil = camberline +- thickness
if p == 0 || H == 0
    xu =  xc;
    zu =  thdis;
    xl =  xc;
    zl = -thdis;
else
    tht(I1,1) = atan((2*H/p)*(-xc(I1)/(c*p)+1));
    tht(I2,1) = atan((2*H/(1-p^2))*(p-(xc(I2)/c)));
    xu = xc         - thdis.*sin(tht);
    zu = camberline + thdis.*cos(tht);
    xl = xc         + thdis.*sin(tht);
    zl = camberline - thdis.*cos(tht);
end