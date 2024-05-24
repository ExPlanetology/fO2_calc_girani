function OutputArray = calculate_IW(T,P)
Conditions=[T,P];
OutputArray(:,1:2)=Conditions(:,1:2);
N = height(Conditions);
for Pair=1:N
    
    T=Conditions(Pair,1);
    P=Conditions(Pair,2);
    
    %Output current temperature and pressure to screen so the user knows the calculation is
    %progressing:
    
    %ScreenOutput = ['Temperature: ' num2str(T),' K; Pressure: ' num2str(P),' GPa'];
    %disp(ScreenOutput);
    
    R=8.314;
    
    
    
    %Calculate G0(T,P) for all components
    % Wustite components and O2 are calculated here.  Those for Fe (the most
    % stable polymorph) are added below.
    
    %G0(i)= a(i)+b(i)*T+c(i)*T*ln(T)+d(i)*T^2+ e(i)*T^2+f(i)*ln(T)+g(i)*T^3
    % i=1 FeO; i=2 FeO1.5; i=3 O2;  i=4 Fe*;
    % G0(4) is added separately from function StableIronPhase, giving G for
    % most stable polymorph of Fe at t and P.
    
    G0=zeros(4,1);
    
    %These parameters are from Hidayat et al. 2015 for FeO and FeO1.5 and
    %Dinsdale (1991) for O2
    
    %Order: FeO, FeO1.5, O2
    a = [-285203.5 -523138 -13137.52];
    b = [274.2455 73.37019 25.32003];
    c= [-49.19444 -26.96809 -33.627];
    d = [-0.004678477 -0.008836071 -0.00119159];
    e = [297568.8	1498519 525809.556];
    f = [574.4469	25471.09 0];
    g = [0           0 0.00000001356111];
    
    if(T<1000)
        
        a(3)=-6961.7445;
        b(3)=-51.0057;
        c(3)=-22.271;
        d(3)=-0.0101977;
        e(3)=-7629.7484;
        f(3)=0;
        g(3)=0.0000000132369;
    end
    %Only the first 3 values of G0 are filled here.
    for i=1:3
        G0(i)= a(i)+b(i)*T+c(i)*T*log(T)+d(i)*T^2+ e(i)/T+f(i)*log(T)+g(i)*T^3;
    end
    
    %  EOS is from Komabayashi, 2014, with V0 of FeO1.5 taken by extrapolation
    %  of wustite stoichiometry-volume data compiled for this project (16.372±0.070)
    %  This is very similar to value from McCammon and Liu 1984 trend (16.425).
    
    EOS=zeros(2,5);
    %             V0,  K0,    K', alpha0,delta0
    EOS(1,:) =  [12.256 149    3.83 4.5e-05 4.25];
    EOS(2,:) =  [16.372 149    3.83 4.5e-05 4.25];
    
    %Mixing parameters for FeO-FeO1.5 from Hidayat et al. 2015
    %Note thaht XFeO=X1; XFeO1.5=X2
    %RT ln gamma(X1) =((q00+2*q10*(1-X2))*X2^2);
    %RT ln gamma(X2) =((1-X2)^2*(q00+q10-2*q10*X2));
    
    q00=-5.94E+04;
    q10=4.27E+04;
    
    % The "j" loop does equation of state calculations for each phase.
    GPressureAdjust=zeros(2,1);
    for j=1:2
        
        PArray=zeros([100,1]);
        VArray=zeros([100,1]);
        
        Pmax=P;
        Peos=0.0001;
        
        % VinetVolumeFunction2 calculates V(T,P) implicitly, calculated fromm
        % the EOS parameters of Komabayashi (2014) JGR.
        
        V=VinetVolumeFunc2(Peos,T,EOS(j,1),EOS(j,2),EOS(j,3),EOS(j,4),EOS(j,5));
        
        PArray(1)=Peos;
        VArray(1)=V;
        
        % the "k" loop fills the V-P  array for numerical integration of VdP
        % It gets looped for each phase, "j"
        
        for k=2:100
            
            Peos=Pmax*(k-1)/(100-1);
            
            V=VinetVolumeFunc2(Peos,T,EOS(j,1),EOS(j,2),EOS(j,3),EOS(j,4),EOS(j,5));
            PArray(k)=Peos;
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
    
    for l=1:2
        
        G0(l)=G0(l)+GPressureAdjust(l);
        
    end
    
    % Now add in G of Fe metal G(4)
    [GFe0, stable]=StableIronPhase(T,P);
    %GFe0 is an array with 5 values, in order
    %fcc, bcc-a, HCP, bcc-d, liquid
    %"stable" gives index of which phase is most stable
    %To supress liquid for IW calculation, the following code is added.
    %Comment out if liquid Fe is desired.
    
    if (stable==5)
        stable=1;
        for x=2:4
            if GFe0(x)<GFe0(stable)
                stable=x;
            end
        end
    end
    
    
    %Now assign stable Fe phase G to G0(4)
    G0(4) = GFe0(stable);
    
    % Now construct delta G of reaction and pass to function that calculates
    % stoichiomestry of wustite in equilibrium with Fe.
    %  2 FeO1.5 + Fe = 3FeO
    % Delta G= 3 G FeO-G Fe - 2 G FeO1.5
    % Added 4/30/2021 - halve reaction stoichiometry to
    % FeO1.5 + 1/2 Fe= 1.5 FeO
    %This should get rid if X2^2 term in expression of log(X1/X2), which seems
    %to allow for negative X2 Solutions.
    
    
    DeltaGFeSat= (3*G0(1)-2*G0(2)-G0(4));
    
    %Equilibrium stoichiometry is given by zero to reaction
    %Note thaht XFeO=X1; XFeO1.5=X2
    % Delta GFeSat+RTln(X1^3/X2^2) + 3 RT ln gamma(X1)-2 RT ln gamma (X2)
    
    %RT ln gamma(X1) =((q00+2*q10*(1-X2))*X2^2);
    %RT ln gamma(X2) =((1-X2)^2*(q00+q10-2*q10*X2));
    
    %Substituting all into one function written as f(X2), the Function that should be zeroed is
    %DeltaGFeSat + R*T *log((1-X2)^3/X2^2) + ...
    %3*((q00+2*q10*(1-X2))*X2^2)- 2*((1-X2)^2*(q00+q10-2*q10*X2)
    %Square it and seek zero with a Newton-Raphson scheme.
    
    func = @(x2)(0.5*DeltaGFeSat+R*T*log((1-x2)^1.5/x2)...
        +0.5*3*((q00+2*q10*(1-x2))*x2^2)-0.5*2*((1-x2)^2*(q00+q10-2*q10*x2)))^2;
    
    %Convergence criterion and initial guess for stoichiometry.
    epsilon1=1e-9;
    x0=0.05;
    
    %NewtRaphZero is internal function specified at bottom of this script
    % It returns equilbirum value of x2 - the mole fraction of XFeO1.5
    
    x2= NewtRaphZero(func,x0,epsilon1);
    
    %Convert to y [Fe(1-y)O]
    
    y=((1-2/(x2+2)));
    x1=1-x2;
    
    %log fO2 is given by reaction FeO + (1/4) O2 = FeO1.5
    % Delta G = GFeO1.5-GFeO-1/4 G O2
    % Delta G = - RT ln K = -RT ln (aFeO1.5/aFeOfO2^1/4)
    % Delta G + RT ln (XFeO1.5/XFeO + RT ln gamma XFeO1.5-RT ln gamma FeO
    
    %X1 is FeO, X2 is FeO1.5
    %gamma1=RT ln gamma(X1) =((q00+2*q10*(1-X2))*X2^2);
    %gamma2=RT ln gamma(X2) =((1-X2)^2*(q00+q10-2*q10*X2));
    
    DeltaGfO2=G0(2)-G0(1)-0.25*G0(3);
    
    gamma1=((q00+2*q10*(1-x2))*x2^2);
    gamma2=((1-x2)^2*(q00+q10-2*q10*x2));
    
    logfO2 =real((4/(2.303*T*8.314)) * (DeltaGfO2 + R*T*log(x2/x1)+ gamma2-gamma1));
    
    OutputArray(Pair,3)=y;
    OutputArray(Pair,4)=logfO2;
    
    
    
end
end
