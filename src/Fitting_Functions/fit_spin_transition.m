% Read the  LS amount at any pressure from Mossbauer data at HP and ambient T
% (PM10-05, BG-08, BG-01) 
clc, clear, close all

in_directory_data  = pwd();                                   %enter the directory where the fitting functions are placed
out_directory_data = fullfile(pwd(), '..');                   %enter the directory where the fo2_calc script is placed
filename_data_LS   = "import_LS_amount.xlsx";
filepath_data_LS   = fullfile(in_directory_data,filename_data_LS);

%% My Data 
Spintransition    = readtable(filepath_data_LS);

%PM10-05
x_PM05            = Spintransition.Pressure_PM05(~isnan(Spintransition.Pressure_PM05));
y_PM05            = Spintransition.LS_Amount_PM05(~isnan(Spintransition.LS_Amount_PM05));

%BG-08
x_BG08            = Spintransition.Pressure_BG08(~isnan(Spintransition.Pressure_BG08));
y_BG08            = Spintransition.LS_Amount_BG08(~isnan(Spintransition.LS_Amount_BG08));

%BG-01
x_BG01            = Spintransition.Pressure_BG01(~isnan(Spintransition.Pressure_BG01));
y_BG01            = Spintransition.LS_Amount_BG01(~isnan(Spintransition.LS_Amount_BG01));

%% Fitting
%PM10-05
[Fit_PM05, gof]   = fit(x_PM05,y_PM05, 'gompertz','Normalize','on')

%BG-08
[Fit_BG08, gof]   = fit(x_BG08,y_BG08, 'gompertz','Normalize','on')

%BG-01
[Fit_BG01, gof]   = fit(x_BG01,y_BG01, 'poly2','Normalize','on')

%% Save fitting
filepath = fullfile(out_directory_data,"fitting_nFeLS.mat");
save(filepath, "Fit_PM05", "Fit_BG08","Fit_BG01")

%% Plot
% PM10-05
figure;
scatter (x_PM05,y_PM05,'MarkerEdgeColor','blue','MarkerFaceColor','blue') 
hold on 

% BG-08
scatter (x_BG08,y_BG08, 'Marker', 's','MarkerEdgeColor','blue','MarkerFaceColor', 'blue') 
hold on 

% BG-01
scatter (x_BG01,y_BG01, 'Marker', 's','MarkerEdgeColor','r','MarkerFaceColor', 'red') 
hold on 

%Figure 1
axis([0 180 -0.05 1.05])
legend('nFe2+ PM10-05','nFe2+ BG-08', 'nFe3+ BG-01');
title('Data from Mossbauer Spectroscopy, This study');
legend('Location','NorthWest');
xlabel('Pressure (GPa)','FontSize',12);
ylabel('nFe LS','FontSize',12);
box on;
% saveas(gcf, fullfile(out_directory_data, 'Data_all.png'));


%% Plot Fit
%PM10-05
figure;
scatter (x_PM05,y_PM05,'Marker', 'o','MarkerEdgeColor','blue','MarkerFaceColor','blue') 
hold on
plot(Fit_PM05,'b')
hold on

% BG-08
scatter (x_BG08,y_BG08,'Marker', 's','MarkerEdgeColor','blue','MarkerFaceColor','blue') 
hold on
plot(Fit_BG08,'b')
hold on
%plot (Fit_BG08,'b',x_BG08,y_BG08, 'bs') %bs=blue 


% BG-01
scatter (x_BG01,y_BG01,'Marker', 'o','MarkerEdgeColor','r','MarkerFaceColor','r') 
hold on
plot(Fit_BG01,'r')
hold on
% %plot(Fit_BG01,'r', x_BG01,y_BG01, 'rs') 

% hold on 
% plot(x_BG01,Fit_BG01,"Color","#EDB120")

%Plot setup
%Figure 2
axis([0 180 -0.05 1.05])
legend ('PM10-05', 'gompertz PM10-05', 'BG-08','gompertz BG-08 ', 'BG-01','Poly2 BG-01',' sigmoidal', ' sigmoidal');
title('Linear Regression Relation Between Fe2+ in the LS state and Pressure');
legend('Location','NorthWest');
xlabel('Pressure (GPa)','FontSize',12);
ylabel('nFe LS','FontSize',12);
box on;

% % saveas(gcf, fullfile(out_directory_data, 'Fitting_LS_amount_all.png'));
% out_filename    = strcat("Fitting_LS","_", string(datetime("today","Format","yyyy_MM_dd")));
% out_file        = fullfile(out_directory_data,out_filename);
% exportgraphics(gcf,strcat(out_file, ".eps"),'BackgroundColor','none','ContentType','vector')
% exportgraphics(gcf,strcat(out_file, ".jpg"),'Resolution',600)





