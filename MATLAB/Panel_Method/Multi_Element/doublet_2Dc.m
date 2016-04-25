function [u,w] = doublet_2Dc(mu,x,z,xL,zL,xR,zR)
%
% [u,w] = doublet_2Dc(mu,x,z,xL,zL,xR,zR)
%
% Calculate the velocity (u,w) at a point (x,z) due to a doublet panel with constant
% strength mu and end points (xL,zL) and (xR,zR).
%

% shift and rotate (xR,zR) and (x,z) so that (xL,zL) is at the origin
% and (xR,zR) is on the x-axis of local panel coordinate system
al_i = atan2(zL-zR,xR-xL);		% angle of panel inclination
xzR = [cos(al_i) -sin(al_i); sin(al_i) cos(al_i)]*[xR-xL zR-zL]';	% (xR,zR) in local panel coordinates
xR = xzR(1); zR = xzR(2);
xz = [cos(al_i) -sin(al_i); sin(al_i) cos(al_i)]*[x-xL z-zL]';		% (x,z) in local panel coordinates
x = xz(1); z = xz(2);

% compute induced velocity components in local panel coordinates
if isinf(xR)
	u = mu/(2*pi)*z/(x^2+z^2);
	w = -mu/(2*pi)*x/(x^2+z^2);
else
	u = mu/(2*pi)*(z/(x^2+z^2)-z/((x-xR)^2+z^2));
	w = -mu/(2*pi)*(x/(x^2+z^2)-(x-xR)/((x-xR)^2+z^2));
end

% rotate velocity components back to airfoil coordinate system
uw = [cos(al_i) sin(al_i); -sin(al_i) cos(al_i)]*[u w]';
u = uw(1); w = uw(2);
end