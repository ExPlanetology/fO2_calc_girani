%Script to calculate a geotherm along an adiabatic thermal profile based on
%previous work of Deng et. al., (2020) and Stixrude et.al., (2009)

clc, clear, close all

in_directory_data           = pwd();                        %enter the directory where the fitting functions are placed
out_directory_data          = fullfile(pwd(), '..');                        %enter the directory where the fo2_calc script is placed
filename_Tprof              = "geotherms";
ext_data_Tprof              = ".xlsx";
filepath_data_Tprof         = fullfile(in_directory_data,strcat(filename_Tprof,ext_data_Tprof));

%% Data 
%Cold and hot Thermal profile from Deng et. al., (2020)
Tprof  = readtable(filepath_data_Tprof);

%Cold Thermal Profile (From 2100K)
x_Tprof_cold       = Tprof.Pressure_2100(~isnan(Tprof.Pressure_2100));
y_Tprof_cold       = Tprof.Temperature_2100(~isnan(Tprof.Temperature_2100));

%Hot Thermal Profile (From 2500K)
x_Tprof_hot       = Tprof.Pressure_2500(~isnan(Tprof.Pressure_2500));
y_Tprof_hot       = Tprof.Temperature_2500(~isnan(Tprof.Temperature_2500));

%% Fitting
%Cold Thermal Profile (From 2100K)
[Fit_Tprof_cold, gof]   = fit(x_Tprof_cold,y_Tprof_cold, 'spline','Normalize','on');

%Hot Thermal Profile (From 2500K)
[Fit_Tprof_hot, gof]   = fit(x_Tprof_hot,y_Tprof_hot, 'spline','Normalize','on');

%% Save fitting
filepath = fullfile(out_directory_data,"geotherm_cold.mat");
save(filepath, "Fit_Tprof_cold")
filepath = fullfile(out_directory_data,"geotherm_hot.mat");
save(filepath, "Fit_Tprof_hot")

%% Plot Fit
%Cold Thermal Profile (2100K)
figure;
scatter (x_Tprof_cold,y_Tprof_cold,'Marker', 'o','MarkerEdgeColor','blue','MarkerFaceColor','blue') 
hold on
plot(Fit_Tprof_cold,'b')
hold on
scatter (x_Tprof_hot,y_Tprof_hot,'Marker', 'o','MarkerEdgeColor','red','MarkerFaceColor','red') 
hold on
plot(Fit_Tprof_hot,'r')
hold on

% Plot settings
legend('Cold','spline', 'Hot','spline');
title('Adibatic temperature profile');
legend('Location','NorthWest');
xlabel('Pressure (GPa)','FontSize',12);
ylabel('Temperature (K)','FontSize',12);

% saveas(gcf, fullfile(out_directory_data, 'Thermal Profile of MO.png'));
% out_filename    = strcat("Geotherms","_", string(datetime("today","Format","yyyy_MM_dd")));
% out_file        = fullfile(out_directory_data,out_filename);
% exportgraphics(gcf,strcat(out_file, ".eps"),'BackgroundColor','none','ContentType','vector')
% exportgraphics(gcf,strcat(out_file, ".jpg"),'Resolution',600)
