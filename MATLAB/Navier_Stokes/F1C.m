function [F1C] = F1C(ue,uw,us,un,vs,vn,dx,dy)
F1C=-((ue.^2)-(uw.^2))/dx-((un.*vn)-(us.*vs))/dy;
end

