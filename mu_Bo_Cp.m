function [ formation_Vol_Factor,Viscosity,Cp_v] = mu_Bo_Cp(T,P)
% Viscosity calculates reservoir oil viscosity
%global Rs gamma_g API gamma_o
Rs=30; gamma_g=.65; API=20; gamma_o=.934;
tem=T-460;
% bubble point using Standing correlation
a=0.00091*tem-.0125*API;
Pb=18.2*(((Rs/gamma_g)^.83).*(10.^a)-1.4);
% Glaso’s Correlation for dead oil viscosity
a=10.313*(log10(tem))-36.447;
mu_dead=(3.141E+10)*(tem.^(-3.444)).*((log10(API)).^a);  % dead oil viscosity
% The Beggs-Robinson Correlation  for saturated oil viscosity
A=10.715*((Rs+100).^(-0.515));
b=5.44*((Rs+150).^(-.338));
mu_bubble=A.*(mu_dead.^b);  % saturated oil viscosity
% The Vasquez-Beggs Correlation for undersaturated oil viscosity
S=(-3.9E-5).*P-5;
m=2.6*(P.^1.187).*(10.^S);
Viscosity=mu_bubble.*((P./Pb).^m); % reservoir oil viscosity
%Vasquez-Beggs’s for oil formation volume factor
Bob=Rs*((gamma_g/gamma_o)^.526)+.968*tem;
A=-6.58511+2.91329*log10(Bob)-.27683*((log10(Bob)).^2);
formation_Vol_Factor_Pb=1+(10.^A);% Glaso's Correlation for oil formation volume factor at bubble point
A=(1E-5)*(-1433+5*Rs+17.2*tem-1180*gamma_g+12.61*API);
formation_Vol_Factor=formation_Vol_Factor_Pb.*exp(-A.*log(P./Pb));
% specific heat 
Cp_v=API*((-1.39E-6).*tem+1.847E-3)+tem.*(6.312E-4)+.352;
end

