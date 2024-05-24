function [G0] = DeltaG_zero(a,b,c,d,e,f,Tf)
% deltaG0 calculation according to Deng et al. 2020
% Return rang of T along which calculate fO2

G0 = a+ b*Tf + c*Tf*log(Tf) + d*Tf^2+ e*Tf^-1 + f*Tf^0.5;

end