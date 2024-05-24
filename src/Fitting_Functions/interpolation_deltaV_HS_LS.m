% Read the (VHS_Fe3-VLS_Fe3) and (VHS_Fe2-VLS_Fe2) from Deng et. al.,
% (2020), (DOI:10.1038/s41467-020-15757-0) and calculate the intermediate value

clc, clear, close all

in_directory_data  = pwd();                                   %enter the directory where the fitting functions are placed
out_directory_data = fullfile(pwd(), '..');                   %enter the directory where the fo2_calc script is placed
filename_data      = "import_deltaHS_LS.xlsx";
filepath_data      =  fullfile(in_directory_data,filename_data);

%% Data
deltaV_LS        = readtable(filepath_data);

x_Temperature    = deltaV_LS.Temperature;
y_deltaV_LS_Fe3  = deltaV_LS.Fe3;
y_deltaV_LS_Fe2  = deltaV_LS.Fe2;

%% Interpolation
[Fit_deltaV_Fe3]   = fit(x_Temperature, y_deltaV_LS_Fe3, 'linear','Normalize','on');  % (VFeO1.5_HS-VFeO1.5_LS)

[Fit_deltaV_Fe2]   = fit(x_Temperature, y_deltaV_LS_Fe2, 'linear','Normalize','on');   % (VFeO_HS-VFeO_LS)

%Fitting 
plot(Fit_deltaV_Fe3,'r');
hold on
plot(Fit_deltaV_Fe2,'b');

%% Save fitting
filepath = fullfile(out_directory_data,"interpolation_deltaV_HS_LS.mat");
save(filepath, "Fit_deltaV_Fe3", "Fit_deltaV_Fe2");

%% Plot setup
legend('FeO_1_._5','FeO' );
xlabel('Temperature (K)');
ylabel('ΔV (V_H_S-V_L_S)');

% % Save plot
% out_filename     = strcat(filename_data,"_", string(datetime("today","Format","yyyy_MM_dd")));
% out_file         = fullfile(out_directory_data,out_filename);
% exportgraphics(gcf,strcat(out_file, ".jpg"),'Resolution',600)

% out_filename_eps = strcat(filename_data, "_", string(datetime("today", "Format", "yyyy_MM_dd")), ".eps");
% out_file_eps     = fullfile(out_directory_data, out_filename_eps);
% print(gcf, out_file_eps, '-depsc', '-r600');



