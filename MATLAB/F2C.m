function [F2C] = F2C(vn,vs,ve,vw,ue,uw,dx,dy)
F2C=-((ue.*ve)-(uw.*vw))/dx-((vn.^2)-(vs.^2))/dy;
end

