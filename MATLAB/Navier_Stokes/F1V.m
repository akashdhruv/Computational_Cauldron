function [F1V] = F1V(uP,uE,uW,uN,uS,dx,dy,mu)
F1V=(mu/dx)*(((uE-uP)/dx)-((uP-uW)/dx))+(mu/dy)*(((uN-uP)/dy)-((uP-uS)/dy));
end

