function [FV] = FV(uP,uE,uW,uN,uS,dx,dy,mu)
FV=(mu/dx)*(((uE-uP)./dx)-((uP-uW)/dx))+(mu/dy)*(((uN-uP)/dy)-((uP-uS)/dy));
end

