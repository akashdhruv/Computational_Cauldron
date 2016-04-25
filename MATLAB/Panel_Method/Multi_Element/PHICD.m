function [phi] = PHICD(mu,x,z,xL,zL,xR,zR)
%
% Calculate the potential phi at a point (x,z) due to a doublet panel with constant
% strength mu and end points (xL,zL) and (xR,zR).
% shift and rotate (xR,zR) and (x,z) so that (xL,zL) is at the origin
% and (xR,zR) is on the x-axis of local panel coordinate system
al_i = atan2(zL-zR,xR-xL);
ala_i= al_i;
xzR = [cos(ala_i) -sin(ala_i); sin(ala_i) cos(ala_i)]*[xR-xL zR-zL]';	% (xR,zR) in local panel coordinates
xR = xzR(1); 
zR = xzR(2);
xz = [cos(ala_i) -sin(ala_i); sin(ala_i) cos(ala_i)]*[x-xL z-zL]';		% (x,z) in local panel coordinates
x = xz(1); 
z = xz(2);
% compute potential in local panel coordinates
phi=-(mu/(2*pi))*((atan2(z,x-xR))-(atan2(z,x)));
end