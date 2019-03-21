function [phi,Der_phi] = phi_Der_phi(P,phi_ref)
%phi_Der_phi calculates porosity and derivative of porosity 
Cf=1E-6;% formation comressiblity
phi=phi_ref.*exp((P-14.7)*Cf);
phi(phi>0.8)=1;
Der_phi=phi*(Cf);
Der_phi(Der_phi==(Cf))=0;
end

