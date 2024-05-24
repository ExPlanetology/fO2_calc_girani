function [V] = VinetVolumeFunc2(P,T,V0,K0,Kprime,alpha0,delta0)

%This function calculates phase volume at P (GPa) and T (K)
%from Vinet-Anderson-Gruneisen EOS of Komabayashi (2014)
%This function is called by script IWCalcHirschmann2.m
%  MMH January 2021
% Modifed, August 2021

%These variables below are now passed by function call
%V0 = 6.82;  fcc iron
%K0 =163.4;   fcc iron
%Kprime = 5.38;  fcc iron
%delta0=5.5;  fcc iron
%alpha0=7e-005;  fcc iron
%T=2000;
%P=50;

T0=298.15;

%This is written as a nested function  so that only one variable (P) is
%passed to  it.  This allows it to be called by other functions such as
%integral.

V=VinetVolumeSubFunc(P);
%Solve for V at P and 298 K (VT0)
% x, VT0,P298 are variables that are only used internally
%Same is true for kappa, used below in thermal expansion calc.

function [V]=VinetVolumeSubFunc(P)


vinet= @(x,P298) 3*K0*(x)^-2*(1-x)*exp(1.5*(Kprime-1)*(1-x))-P298;

% x= (VT0/V0)^(1/3)
% Use fsolve to determine x implicitly
% 0.95 is standard initial guess for x

options = optimset('Display','off');
xsol=fsolve(@(x) vinet(x,P),0.95,options);

VT0_P=V0*xsol^3;

% Then use VT0_P to determine alpha
%Equation 5 from Komabyashi (2014)
kappa=1.4;

alpha=alpha0*exp((-delta0/kappa)*(1-(VT0_P/V0)^kappa));

% Now determine V(T,P) from VT0_P by thermal expansion
V=VT0_P*exp(alpha*(T-T0));
    end

    
end


