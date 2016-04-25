function [n,t] = n_t_vectors(xL,zL,xR,zR)
%
% [n,t] = n_t_vectors(xL,zL,xR,zR)
%
% Compute normal vector n and tangent vector t for a panel with end points (xL,zL) and (xR,zR).

al_i = atan2(zL-zR,xR-xL);			% angle of panel inclination
n = [sin(al_i) cos(al_i)]';			% normal vector
t = [cos(al_i) -sin(al_i)]';		% tangent vector
end