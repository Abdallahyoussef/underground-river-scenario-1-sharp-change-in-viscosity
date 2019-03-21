function [der_formation_Vol_Factor_P] = Der_Bo(T,P)
% Calculating derivative of oil proprties with respect to pressure and temp
global Rs gamma_g API gamma_o
tem=T-460;
A=Rs/gamma_g;
B=(P./((91*10.^((91*tem)./100000 - API./80).*A.^(83/100))./5 - 637./25));
C=log(10);
D=gamma_g/gamma_o;
% derivative of oil formation volume factor with respect to pressure 
der_formation_Vol_Factor_P=-(exp(-log(B).*((1261*API)./10000000 + Rs./20000 - (59.*gamma_g)./5000 + (43.*tem)./250000 - 1433./100000)).*(10.^((820018234901387.*log((121.*tem)./125 + Rs.*(D).^(263./500)))./(281474976710656.*C) - (27683.*log((121.*tem)./125 + Rs.*(D).^(263/500)).^2)./(100000.*C.^2) - 463385920971777./70368744177664) + 1).*((1261.*API)./10000000 + Rs./20000 - (59.*gamma_g)./5000 + (43.*tem)./250000 - 1433./100000))./P;
end

