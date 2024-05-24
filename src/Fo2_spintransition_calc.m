% Calculate the fo2 as a function of P-T and spin transition
% (deltaGf = deltaG0 + integral(deltaV*dP)).
%
% MO composition from McDonough and Sun, 1995
% (https://doi.org/10.1016/0009-2541(94)00140)
%
% Fe3+/Fetot from 0.037 and 0.1 (see Hirschmann 2022, DOI:10.1016/j.gca.2022.04.005).
%
% Activity calculated assuming ln_gamma equal to 0 (Sossi et.al., 2020 et.al., 2020, DOI: 10.1126/sciadv.abd1387)
% or based on Margules parameters from Deng et. al., (2020) (DOI:10.1038/s41467-020-15757-0). 
% The user can choose one or the other option using a flag in the code.
% 
% The user can choose to calculate the fo2 at a specific temperature (0) or
% along a getherm (1) using the geotherm shown by Deng et. al., (2020)
% 
% Delta G0 from Deng et. al., (2020) (equation in supplementary note 1).
%
% Thermal pressure EoS (basd on a BM4) is implemented as in eq. 2 of Deng et. al., (2020).
% The parameters used in the EoS are from Hirschmann (2022).
%
% The nFe2+ and nFe3+ are obtained from Synchrotron Mossbauer data on PM10-05,
% BG-08, BG-01 obtained at 298K. The dependece of nFe2+ and nFe3+ with
% respect to P [GPa] is fitted in another script. The fitting functions are
% then loaded in this script.
%
% The nFe2+ and nFe3+ are obtained from measurements at 298K but they are
% used for calculations at various P,T in this script. In order to
% determine the appropriate nFe2+ and nFe3+ at each P,T
% we do the following steps (example only for Fe3+, Fe2+ is analogous):
% 1) the VFe3+ at the given P,T is calculated from the
% corresponding EoS, i.e VFe3+ = EoS(P,T)
% 2) we use the EoS to calculate the Pscale = EoS(VFe3+, 298K)
% where VFe3+ is the volume obtained at the previous step.
% 3) The nFe is obtained from the fitting function at the Pscale obtained
% in the previous step.
%
% This approach assumes that the value of nFe2+ and nFe3+ only depends on volume.


clc, clear, close all
%% Input Parameters
directory_data      = pwd();                     %enter the directory where thw fo2_calc script is placed
out_directory_data  = pwd();                     %enter the desired directory for the output file
filename_data       = "Input_Parameters_All.xlsx";
nFe_fitting_file    = "fitting_nFeLS.mat";
geotherm            = "geotherm_cold.mat";                  %choose file name according to the geotherm
deltaV_HS_LS        = "interpolation_deltaV_HS_LS.mat";

print               = true;                                 %If true, print the results to the console


N_Avogadro	        = 6.022e+23;
R                   = 8.314;        % [J/mol*K] Gas costant
Fe3_Fetot           = 0.065;        % choose the Fe3+/Fetot
Tflag               = 1;            % choose (0) to calculate only at a specific temperature or (1) to calculate along the geotherm
Tchoice             = 2432;         % [K] Set this temperature when Tflag==0. The calculation will be performed at this temperature.
Pmin                = 0;            % [GPa] minimum P along the geotherm
Pmax                = 140;          % [GPa] maximum P along the geotherm
nP                  = 10;           % choose number of steps along the geotherm
lngamma_Flag        = 0;            % If (0) From Sossi et.al., (2020), if (1) calculate value from Deng et.al., (2020)'s margules parameters
P0                  = 0;            %[GPa] Reference pressure for thermodynamic calculation(i.e. reference P for the calculation of volume integrals and for the Gibbs energy)


% %12.5 mol% Fe Deng et. al., (2020)
%Eos Fe2+
% params_Fe2 = struct;
% params_Fe2.V0  = 1180.10;         %[A3]
% params_Fe2.K   = 26.76;           %[GPa]
% params_Fe2.Kp  = 2.8;
% params_Fe2.Kpp = 0.01;            %[GPa-1]
% params_Fe2.a   = 35.70;
% params_Fe2.b   = 71.10;
% params_Fe2.c   = 36.60;
% params_Fe2.T0  = 3000;            %[K]

% Eos Fe3+
% params_Fe3.V0  = 1204.69;         %[A3]
% params_Fe3.K   = 23.18;           %[GPa]
% params_Fe3.Kp  = 3.22;
% params_Fe3.Kpp = 0.01;            %[GPa-1]
% params_Fe3.a   = 34.53;
% params_Fe3.b   = 68.64;
% params_Fe3.c   = 35.27;
% params_Fe3.T0  = 3000;            %[K]

%12.5 mol% Fe (Hirschmann, 2022)
% Eos Fe2+
params_Fe2 = struct;
params_Fe2.V0  = 1180.114;          %[A3]
params_Fe2.K   = 26.947;            %[GPa]
params_Fe2.Kp  = 2.803;
params_Fe2.Kpp = 0.012;
params_Fe2.a   = 35.794;
params_Fe2.b   = 71.103;
params_Fe2.c   = 36.595;
params_Fe2.T0  = 3000;              %[K]

% Eos Fe3+
params_Fe3.V0  = 1204.764;          %[A3]
params_Fe3.K   = 23.195;            %[GPa]
params_Fe3.Kp  = 3.216;
params_Fe3.Kpp = 0.009;             %[GPa-1]
params_Fe3.a   = 34.526;
params_Fe3.b   = 68.644;
params_Fe3.c   = 35.271;
params_Fe3.T0  = 3000;              %[K]

% Parameters for deltaGzero calculation

deltag0.a   = -3.310*10^5;     %From Deng et. al., (2020), supplementary note 1
deltag0.b   = -190.379;
deltag0.c   = 14.785;
deltag0.d   = -1.649*10^-3;
deltag0.e   = 9.348*10^6;
deltag0.f   = 1.077*10^4;

%% Initialize
% Make a table to save all the results of the calculation
Pgeotherm = linspace(Pmin,Pmax,nP);
PTdata    = table(zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),...
    zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),...
    zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),...
    zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),zeros(size(Pgeotherm,2),1),...
   'VariableNames', ...
    {'P', 'T', 'Fe3_Fetot','deltaG0','lngamma','Fe3_HS_LS', 'Fe2_HS_LS','nFe3','nFe2_PM05','nFe2_BG08',...
    'DeltaGf_HS','DeltaGf_PM10_05','DeltaGf_BG_08','IW',...
    'log_fO2_HS','log_fO2_PM10_05','log_fO2_BG_08',...
    'deltaIW_HS','deltaIW_PM10_05', 'deltaIW_BG_08'});

%% Import files
filepath = fullfile(directory_data,filename_data);
data = importfile_fo2(filepath, 1, [2,11]);

%Import fitting of nFeLS
filepath = fullfile(directory_data,nFe_fitting_file);
load(filepath);

%Import fitting of geotherms
filepath        = fullfile(directory_data,geotherm);
temp            = struct2cell(load(filepath));      
Fit_Tgeotherm   = temp{1};

%import deltaV_LS for FeO1.5 (i.e.(VFe3_HS_VFe3_LS) and FeO (i.e.(VFe2_HS_VFe2_LS)
filepath =  fullfile(directory_data,deltaV_HS_LS);
load(filepath)

%% Calculte Moles for each oxide

Moles_Si            = data.Weight(1)/(data.AtomicMass(1)+(data.AtomicMass(10)*2));
Moles_Al            = data.Weight(2)/((data.AtomicMass(2)*2)+(data.AtomicMass(10)*3));
Moles_Fe            = data.Weight(3)/(data.AtomicMass(3)+data.AtomicMass(10));
Moles_Mg            = data.Weight(4)/(data.AtomicMass(4)+data.AtomicMass(10));
Moles_Ca            = data.Weight(5)/(data.AtomicMass(5)+data.AtomicMass(10));
Moles_K             = data.Weight(6)/((data.AtomicMass(6)*2)+data.AtomicMass(10));
Moles_Na            = data.Weight(7)/((data.AtomicMass(7)*2)+data.AtomicMass(10));
Moles_Ti            = data.Weight(8)/(data.AtomicMass(8)+(data.AtomicMass(10)*2));
Moles_P             = data.Weight(9)/((data.AtomicMass(9)*2)+(data.AtomicMass(10)*5));
Moles_Tot           = Moles_Si + Moles_Al + Moles_Fe + Moles_Mg +  Moles_Ca  + Moles_K +  Moles_Na + Moles_Ti + Moles_P;

%% Calculte Molar Fraction and Fe3+ and Fe2+ contribution
MolarF_Si  = Moles_Si/Moles_Tot;
MolarF_Al  = Moles_Al/Moles_Tot;
MolarF_Fe  = Moles_Fe/Moles_Tot;
MolarF_Mg  = Moles_Mg/Moles_Tot;
MolarF_Ca  = Moles_Ca/Moles_Tot;
MolarF_K   = Moles_K/Moles_Tot;
MolarF_Na  = Moles_Na/Moles_Tot;
MolarF_Ti  = Moles_Ti/Moles_Tot;
MolarF_P   = Moles_P/Moles_Tot;

Fe3_amount     = (MolarF_Fe)*(Fe3_Fetot);
Fe2_amount     = (MolarF_Fe)-(Fe3_amount);

%%

    for iP = 1:length(Pgeotherm)
        Pcurrent = Pgeotherm(iP);
        if Tflag == 0
           disp('T is zero, allocated Tcurrent');
           Tcurrent = Tchoice;
        else 
           disp(['T is not zero, allocated Tcurrent: ', num2str(Fit_Tgeotherm(Pcurrent))]);
           Tcurrent = Fit_Tgeotherm(Pcurrent);
        end
    
        %% Calculte Gibbs free energey at ambient conditions
        deltaG0  = DeltaG_zero(deltag0.a, deltag0.b,deltag0.c,deltag0.d, deltag0.e, deltag0.f, Tcurrent);     %[J/mol] From Deng et. al., (2020), supplementary note 1

        %% Calculate deltaV over the range between P0 and the current P
        Prange_integral         = linspace(P0,Pcurrent,1000);         % [GPa] Discretization of P to calculate the deltaV and the integral
        V_Fe2                   = zeros(size(Prange_integral));
        V_Fe3                   = zeros(size(Prange_integral));
        deltaV_HS               = zeros(size(Prange_integral));
        deltaV_LS_PM10_05       = zeros(size(Prange_integral));
        deltaV_LS_BG_08         = zeros(size(Prange_integral));
        deltaV_PM10_05          = zeros(size(Prange_integral));
        deltaV_BG_08            = zeros(size(Prange_integral));
        Vinit                   = 1180; %[A3]
    
        for i = 1:length(Prange_integral)
    
            %------DeltaV of HS------
            V_Fe2(i)        = BM_Th_VPT(Prange_integral(i),Vinit,Tcurrent,params_Fe2.V0,params_Fe2.K,params_Fe2.Kp,params_Fe2.Kpp,params_Fe2.a,params_Fe2.b,params_Fe2.c,params_Fe2.T0); %[A^3]
            V_Fe3(i)        = BM_Th_VPT(Prange_integral(i),Vinit,Tcurrent,params_Fe3.V0,params_Fe3.K,params_Fe3.Kp,params_Fe3.Kpp,params_Fe3.a,params_Fe3.b,params_Fe3.c,params_Fe3.T0); %[A^3]
            deltaV_HS(i)    = (V_Fe3(i) - V_Fe2(i))/2*N_Avogadro*1e-30; %[m3/mol]
    
            %%PM10-05
            %------DeltaV of spin transition------
            % Find the P at 298K that corresponds to the V at the current P and T
            Vinit                       = 1180; %[A^3]
            V                           = BM_Th_VPT(Prange_integral(i),Vinit,Tcurrent,params_Fe2.V0,params_Fe2.K,params_Fe2.Kp,params_Fe2.Kpp,params_Fe2.a,params_Fe2.b,params_Fe2.c,params_Fe2.T0);      %[A^3]
            Pscaled                     = BM_Th_PVT(V,298,params_Fe2.V0,params_Fe2.K,params_Fe2.Kp,params_Fe2.Kpp,params_Fe2.a,params_Fe2.b,params_Fe2.c,params_Fe2.T0);                   %[GPa]
            %Get the nFe2+ from the fitting
            nFe2_PM10_05                = Fit_PM05(Pscaled);

            % Find the P at 298K that corresponds to the V at the current P and T 
            V                           = BM_Th_VPT(Prange_integral(i),Vinit,Tcurrent,params_Fe3.V0,params_Fe3.K,params_Fe3.Kp,params_Fe3.Kpp,params_Fe3.a,params_Fe3.b,params_Fe3.c,params_Fe3.T0);      %[A^3]
            Pscaled                     = BM_Th_PVT(V,298,params_Fe3.V0,params_Fe3.K,params_Fe3.Kp,params_Fe3.Kpp,params_Fe3.a,params_Fe3.b,params_Fe3.c,params_Fe3.T0);                   %[GPa]
            %Get the nFe3 from the fitting
            nFe3                       = Fit_BG01(Pscaled);


            %Get the V(HS-LS) for Fe3+ and Fe2+ from the fitting (done in
            %another script)
            VHS_VLS_FeO1_5              = Fit_deltaV_Fe3(Tcurrent);      %Import the linear function determined for Fe3+ to determinde the deltaV OF Fe3+ between HS and LS configuration (Data at 3000K and 4000K from Deng et. al., 2020)
            VHS_VLS_FeO                 = Fit_deltaV_Fe2(Tcurrent);      %Import the linear function determined for Fe2+ to determinde the deltaV OF Fe3+ between HS and LS configuration (Data at 3000K and 4000K from Deng et. al., 2020)

            %calculate the deltaV of transition        
            deltaV_LS_PM10_05(i)        = - nFe3*(VHS_VLS_FeO1_5) + nFe2_PM10_05*(VHS_VLS_FeO); %[m3/mol]

            %------total DeltaV------
            deltaV_PM10_05(i)           = deltaV_HS(i) + deltaV_LS_PM10_05(i); %[m3/mol]


            %%BG_08
            %------DeltaV of spin transition------
            % Find the P at 298K that corresponds to the V at the current P and T
            Vinit                       = 1180; %[A^3]
            V                           = BM_Th_VPT(Prange_integral(i),Vinit,Tcurrent,params_Fe2.V0,params_Fe2.K,params_Fe2.Kp,params_Fe2.Kpp,params_Fe2.a,params_Fe2.b,params_Fe2.c,params_Fe2.T0);      %[A^3]
            Pscaled                     = BM_Th_PVT(V,298,params_Fe2.V0,params_Fe2.K,params_Fe2.Kp,params_Fe2.Kpp,params_Fe2.a,params_Fe2.b,params_Fe2.c,params_Fe2.T0);                   %[GPa]
            %Get the nFe2 from the fitting
            nFe2_BG_08                  = Fit_BG08(Pscaled);
    
            % Find the P at 298K that corresponds to the V at the current P and T
            V                           = BM_Th_VPT(Prange_integral(i),Vinit,Tcurrent,params_Fe3.V0,params_Fe3.K,params_Fe3.Kp,params_Fe3.Kpp,params_Fe3.a,params_Fe3.b,params_Fe3.c,params_Fe3.T0);      %[A^3]
            Pscaled                     = BM_Th_PVT(V,298,params_Fe3.V0,params_Fe3.K,params_Fe3.Kp,params_Fe3.Kpp,params_Fe3.a,params_Fe3.b,params_Fe3.c,params_Fe3.T0);                   %[GPa]
            %Get the nFe3 from the fitting
            nFe3                        = Fit_BG01(Pscaled);
    
            %calculate the deltaV of transition
            deltaV_LS_BG_08(i)          = - nFe3*(VHS_VLS_FeO1_5) + nFe2_BG_08*(VHS_VLS_FeO); %[m3/mol]    
    
            %------total DeltaV------
            deltaV_BG_08(i)             = deltaV_HS(i) + deltaV_LS_BG_08(i); %[m3/mol]
    
        end
        %% Integrate deltaV in Pressure
        %convert GPa to Pa
        Prange_Pa = Prange_integral*10^9;
    
        %integrate deltaV in dP
        deltaVHS_dP           = trapz(Prange_Pa,  deltaV_HS);         %[m3/mol * Pa] integral of deltaV between Fe3+ - Fe2+ at high spin
        deltaV_dP_PM10_05     = trapz(Prange_Pa,  deltaV_PM10_05);    %[m3/mol * Pa] integral of total deltaV including spin transition for sample : PM10-05
        deltaV_dP_BG_08       = trapz(Prange_Pa,  deltaV_BG_08);      %[m3/mol * Pa] integral of total deltaV including spin transition for sample : BG-08
       
        %% Calculate Gamma(FeO1.5-FeO), X(FeO1.5-FeO)

        if lngamma_Flag == 0
            disp('lngamma_Flag is zero, allocated 0');
            lngamma = 0;                                                       % From Sossi et.al., (2020)
        elseif lngamma_Flag == 1                                               % From Deng et. al., (2020)
            disp('lngamma_Flag is 1, allocated 1');
            lngamma  =  ((MolarF_Si*data.Margules(2)+MolarF_Al*data.Margules(3)+MolarF_Mg*data.Margules(1)+MolarF_Ca*data.Margules(4)+MolarF_K*data.Margules(6))/(R*Tcurrent))+((Fe2_amount-Fe3_amount)*data.Margules(9))/(R*Tcurrent);
            % From Deng et. al., (2020), supplementary note 2
        end

        X_Fe             = log((Fe3_amount/Fe2_amount));
        gamma            = exp(lngamma);
        activity         = gamma*((Fe3_amount/Fe2_amount));                     % Activity(FeO1.5-FeO)

        %% Delta Gibbs free energy final
        %Assign the variable to the table 
        PTdata.P(iP)                    = Pcurrent;
        PTdata.T(iP)                    = Tcurrent;
        PTdata.Fe3_Fetot(iP)            = Fe3_Fetot;
        PTdata.deltaG0(iP)              = deltaG0;
        PTdata.lngamma(iP)              = lngamma;
        PTdata.Fe3_HS_LS(iP)            = VHS_VLS_FeO1_5;    %This is the volume difference of FeO1.5 and FeO between the High spin and low spin state from previous work: Karki et.al., (2018, DOI:10.1029/2018GL077149) and Deng et. al., (2020) (i.e. (VFe3+HS-Fe3+LS) and (VFe2+HS-Fe2+LS)
        PTdata.Fe2_HS_LS(iP)            = VHS_VLS_FeO;
        PTdata.nFe3(iP)                 = nFe3;
        PTdata.nFe2_PM05(iP)            = nFe2_PM10_05;
        PTdata.nFe2_BG08(iP)            = nFe2_BG_08;
        
        
        % Fugacity accounting only for volume changes due to valence state (VFe3+ HS -VFe2+ HS)
        deltaGfinal_HS              = deltaG0 + deltaVHS_dP;                     %[J/mol]
        ln_k_HS                     = deltaGfinal_HS/(R*Tcurrent);
        ln_fo2_HS                   = 4*(ln_k_HS + lngamma + X_Fe);
        fo2_HS                      = exp(ln_fo2_HS);
        log_fO2_HS                  = log10(fo2_HS);
        % Assign variable to the table
        PTdata.DeltaGf_HS(iP)       = deltaGfinal_HS;
        PTdata.log_fO2_HS(iP)       = log_fO2_HS;
    
        % PM10-05 Fugacity accounting also for volume changes due to spin transition
        deltaGfinal_total_PM10_05           = deltaG0 + deltaV_dP_PM10_05;                       %[J/mol]
        ln_k_total                          = deltaGfinal_total_PM10_05/(R*Tcurrent);
        ln_fo2_total                        = 4*(ln_k_total + lngamma + X_Fe);
        fo2_total_PM10_05                   = exp(ln_fo2_total);
        log_fO2_total_PM10_05               = log10(fo2_total_PM10_05);
        % Assign variable to the table
        PTdata.DeltaGf_PM10_05(iP)          = deltaGfinal_total_PM10_05;
        PTdata.log_fO2_PM10_05(iP)          = log_fO2_total_PM10_05;
        
    
        % BG_08 Fugacity accounting also for volume changes due to spin transition
        deltaGfinal_total_BG_08             = deltaG0 + deltaV_dP_BG_08;                       %[J/mol]
        ln_k_total                          = deltaGfinal_total_BG_08 /(R*Tcurrent);
        ln_fo2_total                        = 4*(ln_k_total + lngamma + X_Fe);
        fo2_total_BG_08                     = exp(ln_fo2_total);
        log_fO2_total_BG_08                 = log10(fo2_total_BG_08 );
        % Assign variable to the table
        PTdata.DeltaGf_BG_08(iP)            = deltaGfinal_total_BG_08;
        PTdata.log_fO2_BG_08(iP)            = log_fO2_total_BG_08;
        

    
        %% Calculate IW using the script written by Hirschmann 2021 (DOI: 10.1016/j.gca.2021.08.039)
    
        OutputArray            = calculate_IW(Tcurrent,Pcurrent);
        IWfo2                  = OutputArray(:,4);            % 1Temperature, 2Pressure, 3 .., 4 logfo2 of IW
        PTdata.IW(iP)          = IWfo2;
        
        %% Calculate ΔIW (logfo2-IW)
    
        deltaIW_HS                              = log_fO2_HS - IWfo2;
        deltaIW_withLS_PM10_05                  = log_fO2_total_PM10_05 - IWfo2;
        deltaIW_withLS_BG_08                    = log_fO2_total_BG_08 - IWfo2;

    
        PTdata.deltaIW_HS(iP)                   = deltaIW_HS;
        PTdata.deltaIW_PM10_05(iP)              = deltaIW_withLS_PM10_05; 
        PTdata.deltaIW_BG_08(iP)                = deltaIW_withLS_BG_08;
       
       
        %% Print
        if print
            fprintf("-------------------------------------------------------------------------------------\n")
            fprintf("Iteration %.0f\n", iP)
            fprintf("------Values of model------\n")
    
            fprintf("P [GPa]: %.3f, T [K]: %.3f\n", Pcurrent, Tcurrent)
            fprintf("Fe3+/Fetot: %.3f\n", (Fe3_amount/(Fe2_amount+Fe3_amount)))
    
            fprintf("------Results------\n")
    
            fprintf("XFe: %.3f\n", (Fe3_amount/Fe2_amount))
            fprintf("log(XFe): %.3f\n", X_Fe)
            fprintf("DeltaG0 [J/mol]: %.3f\n", deltaG0)
            %fprintf("DeltaV at high spin: %.3f\n", deltaV_HS)
            % fprintf("deltaV_LS_PM10_05: %.3f\n", deltaV_LS_PM10_05);
            % fprintf("deltaV_LS_BG_08: %.3f\n", deltaV_LS_BG_08);
            % fprintf("DeltaV total: %.3f\n", deltaV)
            fprintf("Integral of deltaV in dP [m3/mol * Pa]: %.5f\n", deltaVHS_dP)
            fprintf("Integral of deltaV in dP (including spin transition) [m3/mol * Pa]: %.5f \n", deltaV_dP_PM10_05)
            fprintf("DeltaG at P and T [J/mol]: %.3f\n", deltaGfinal_HS)
            fprintf("DeltaG at P and T (including spin transition) [J/mol]: %.3f\n", deltaGfinal_total_PM10_05)
    
            fprintf("------oxygen fugacity------\n")
    
            fprintf("logfO2_HS: %.3f\n", log_fO2_HS)
            fprintf("logfO2_IW_: %.3f\n", IWfo2)
            fprintf("deltaIW_HS: %.3f\n", deltaIW_HS)
    
            fprintf("------oxygen fugacity including spin transition------\n")
    
          
            fprintf("logfO2_PM10_05: %.3f\n", log_fO2_total_PM10_05)
            fprintf("logfO2_BG_08: %.3f\n", log_fO2_total_BG_08)
    
        
            fprintf("deltaIW_PM10_05: %.3f\n", deltaIW_withLS_PM10_05)
            fprintf("deltaIW_BG_08: %.3f\n", deltaIW_withLS_BG_08)
    
            % deltaV_LS_PM10_05 = deltaV_LS_PM10_05(end)*10^6;
            % deltaV_LS_BG_08   = deltaV_LS_BG_08 (end)*10^6;
        end
    end

%% Save table
filename_table = 'fo2_calc.xlsx';
out_file        = fullfile(out_directory_data, filename_table);
writetable(PTdata, filename_table, 'Sheet',1)






