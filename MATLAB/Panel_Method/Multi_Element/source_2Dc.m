function [u,w] = source_2Dc(sig,x,z,xL,zL,xR,zR)
%
% [u,w] = doublet_2Dc(mu,x,z,xL,zL,xR,zR)
%
% Calculate the velocity (u,w) at a point (x,z) due to a source panel with constant
% strength sigma and end points (xL,zL) and (xR,zR).
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
	u = sig/(2*pi)*z/(x^2+z^2);
	w = sig/(2*pi)*x/(x^2+z^2);
else
	u = (sig/(4*pi))*log((x^2+z^2)/((x-xR)^2+z^2));
	w = (sig/(2*pi))*(atan2(z,(x-xR))-atan2(z,x));
end

% rotate velocity components back to airfoil coordinate system
uw = [cos(al_i) sin(al_i); -sin(al_i) cos(al_i)]*[u w]';
u = uw(1); w = uw(2);
end