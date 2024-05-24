function [P] = BM_Th_PVT(V,T,V0,K,Kp,Kpp,a,b,c,T0)
% Birch Murnaghan + thermal pressure according to Deng et al., 2020
% Return P given V and T


P = BM_PV(V,V0,K,Kp,Kpp) + (a-b*(V/V0)+c*(V/V0)^2)/1000 * (T-T0);





end

