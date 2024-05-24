function [V] = BM_Th_VPT(P,Vinit,T,V0,K,Kp,Kpp,a,b,c,T0)
% Birch Murnaghan + thermal pressure according to Deng et al., 2020
% Return V[A3] given P and T


func = @(Vinit) BM_Th_PVT(Vinit,T,V0,K,Kp,Kpp,a,b,c,T0) - P; 
V  = fzero(func,Vinit);

end

