function [G, stable]= StableIronPhase(T,P)
%MMH April, 2021, modified July, 2021

%NOTE: THis function reads data from a file each time it is called.  
% The file is SGTEIron.xlsx, which must be present in an appropriate directory for
%the script to run properlly.
%The script would run faster if these were transfered to internal arrays.

%Function calculates the stable polymorph of iron (bcc, fcc, or hcp) at a given T and P and
%passes back its Gibbs free energy (G) and its identity (IronPhase)

%100 kPa values for phases of iron are from the SGTE database (Dinsdale,
%1991, CALPHAD) and are consistent with values adopted by Sundman (1991)
%and Hidayat et al. (2015).  Note that there are two bcc phases - alpha
%and delta and that the databse includes Fe liquid, which is not used in
%the present script.

%EOS for fcc and hcp from Komabayashi (2014).  bcc assumed to have the same
%EOS as fcc except for a larger V0, estimated here to make phase diagram
%accurate.

if T<1811
TableInput =readtable('SGTEIron.xlsx','Sheet','LowT','ReadVariableNames',false);
else
TableInput =readtable('SGTEIron.xlsx','Sheet','HiT','ReadVariableNames',false);

end


%EOS parameters: Order: fcc bcc-alpha HCP bcc-delta liquid
%EOS from Komabayashi, 2014, with bcc the same as fcc, except 
%that V0 of bcc is taken from Dorogokupets et al. 2017

%             V0,    K0,   K',  alpha0,delta0
EOSFe(1,:) =  [6.82   163.4 5.38 7e-05   5.5];
EOSFe(2,:) =  [7.092   163.4 5.38 7e-05   5.5];
EOSFe(3,:) =  [6.753  163.4 5.38 5.8e-05 5.1];
EOSFe(4,:) =  [7.092   163.4 5.38 7e-05   5.5];
EOSFe(5,:) =  [6.88   148   5.8  9.0e-05 5.1];

S = vartype('numeric');
Gparams= TableInput{:,S};

N=width(Gparams);


a=Gparams(1,:);
b=Gparams(2,:);
c=Gparams(3,:);
d=Gparams(4,:);
e=Gparams(5,:);
f=Gparams(6,:);
g=Gparams(7,:);
h=Gparams(8,:);

%gg loop calculates G for each phase

G=zeros(N);

for gg= 1:N
    G(gg)=a(gg) + b(gg)*T + c(gg)*T*log(T) + d(gg)*T^2 + e(gg)*T^3 + f(gg)/T + g(gg)*T^7 + h(gg)*T^-9;
end
    
%Adjust bcc-alpha for magnetic contribution.  Magnetic contribution for fcc
% is << 1 J at temperatures of interest and is ignored

Tc = 1043;
Pfactor = 0.4;
Beta = 2.22;
A = 1.55828482;

if(T<Tc)
   Gmag=1-(1/A)*((79*(T/Tc)^-1/(140*Pfactor))+(474/497)*((1/Pfactor)-1)*((T/Tc)^3/6+(T/Tc)^9/135+(T/Tc)^15/600));
else    
  Gmag = (-1/A)*((T/Tc)^-5/10+(T/Tc)^-15/315+(T/T)^-25/1500);
end

Gmag=Gmag*(8.314*T*log(Beta+1));

G(2)=G(2)+Gmag;


PArray=zeros([100,1]);
VArray=zeros([100,1]);
GPressureAdjust=zeros(5,1);

Pmax=P;


% The "j" loop does equation of state calculations for each phase.
for j=1:5


P=0.0001;
% VinetVolumeFunction2 calculates V(T,P) implicitly, calculated fromm
% the EOS parameters of Komabayashi (2014) JGR.

V=VinetVolumeFunc2(P,T,EOSFe(j,1),EOSFe(j,2),EOSFe(j,3),EOSFe(j,4),EOSFe(j,5));

PArray(1)=P;
VArray(1)=V;

% the "k" loop fills the V-P  array for numerical integration of VdP
% It gets looped for each phase, "j"

for k=2:100
    
    P=Pmax*(k-1)/(100-1);

V=VinetVolumeFunc2(P,T,EOSFe(j,1),EOSFe(j,2),EOSFe(j,3),EOSFe(j,4),EOSFe(j,5));
PArray(k)=P;
VArray(k)=V;
end

VdP=trapz(PArray,VArray);

%This is the VdP integral, calculated numericaly with the trapz function
% from V(T,0) to V(T,P)
%,in (cm^3/mol )*GPa 
%Units of P = GPa, Units of V cm^3/mole
% 1 cm^3/mole= 0.1 J/bar/mole=1000 J/GPa/mole
GPressureAdjust(j)=VdP*1000;


end

for l=1:5
    
    G(l)=G(l)+GPressureAdjust(l);
    
end

stable=1;
for x=2:5
if G(x)<G(stable)
stable=x;
end
end

end




