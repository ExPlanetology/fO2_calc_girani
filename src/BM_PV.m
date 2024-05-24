function [P] = BM_PV(V,V0,K,Kp,Kpp)
%Birch Murnaghan, return P[GPa] given V. The units of V and V0 must be the
%same.

fe = ((V0/V)^(2/3)-1)/2;
P = 3 * K * fe * ((1+2*fe)^(5/2)) * (1+3/2*(Kp-4)*fe +3/2*(K*Kpp+(Kp-4)*(Kp-3)+35/9)*fe^2);
end

