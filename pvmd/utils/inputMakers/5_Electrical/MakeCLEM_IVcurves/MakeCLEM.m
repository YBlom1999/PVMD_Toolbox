% This script makes a Calibrated Lumped Element Model of the solar cell.
% It requires the measured/simulated IV curve of a cell for different
% temperatures and irradiance levels.
%
% For each IV curve, a fit is made with the one-diode equivalent circuit
% model.
% Then, a fit is made between the value of equivalent circuit parameters
% w.r.t. the varying parameter (i.e. temperature absorbed current).
%
% This script will automatically save the results into the correct data
% format such that it can be used by the Toolbox.
% To use this script, fill in the correct data in the input section.
% Make sure that the filenames contain the same structure as the one in the
% example (i.e. 'filename'-_25deg_1000Wm2-'extension')
%
% Developed by Malte Vogt and Youri Blom

%% Input
%The name of the folder where the IV curves are saved
Folder_name = 'Example';

%The name of the IV curve files
JV_name = 'JV_Example';

%The name of the Absorption files (name it '-' if no absorption
%files are present)
Abs_name = 'Abse_Example';

%The column index of the absorption file that represents the absorber
%material. Wavelength should be first column. (name it nan if no absorption
%files are present)
Abs_ind = 6;

%The extension of all files
Extension = '.dat';

%The different temperatures for which the IV curves are measured/simulated.
T_levels = [-35,-25,-15,-5,5,15,25,35,45,55,65,75,85,95];

%The irradiance level during the different temperature measurements/simulations
Irr_fixed = 1000;

%The different irradiance levels for which the IV curves are measured/simulated.
Irr_levels = [50,100,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500];

%The temperature during the different irradiance measurements/simulations
T_fixed = 25;

%Correction for the voltage or current axis (they both need to be positive)
Flip_voltage = 0;
Flip_current = 1;

%The filename of the resulting matlab file
Filename_save = 'Example.mat';

%% Check if files are read properly
%Combine inputs
Filenames.Folder_name = Folder_name;
Filenames.JV_name = JV_name;
Filenames.Abs_name = Abs_name;
Filenames.Extension = Extension;

[Check_Irr,Voltage_Irr,Current_Irr,Jph_Irr] = Check_Curves(Irr_levels,Filenames,Abs_ind,T_fixed,Flip_voltage,Flip_current);
[Check_T,Voltage_T,Current_T,Jph_T] = Check_Curves(Irr_fixed,Filenames,Abs_ind,T_levels,Flip_voltage,Flip_current);

if ~Check_Irr || ~Check_T
    error('Curves are not correct')
end


%% Calculate dependencies

[J_sc_J,J_0_J,R_s_J,R_sh_J,n_J] = Fit_Irr_curves(Voltage_Irr, Current_Irr,Jph_Irr,Irr_levels,T_fixed);
[J_sc_T,J_0_T,R_s_T,R_sh_T,n_T] = Fit_Temp_curves(Voltage_T, Current_T,T_levels);

save(Filename_save,'J_sc_J','J_0_J','R_s_J','R_sh_J','n_J','J_sc_T','J_0_T','R_s_T','R_sh_T','n_T');

%% Functions
function [Check,Voltage,Current,Jph] = Check_Curves(Irr_levels,Filenames,Abs_ind,T_levels,Flip_voltage,Flip_current)
%This script asks the user to check if the IV curves are read and processed
%correctly.
Folder_name = Filenames.Folder_name;
JV_name = Filenames.JV_name;
Abs_name = Filenames.Abs_name;
Extension = Filenames.Extension;

JV_Filename = fullfile(Folder_name,append(JV_name,'_',num2str(T_levels(1)),'deg_',num2str(Irr_levels(1)),'Wm2',Extension));
data = importdata(JV_Filename);
N_points = length(data.data(:,1));
N_curves = max(length(Irr_levels),length(T_levels));
if N_curves == length(Irr_levels)
    Type = 'Irr';
else
    Type = 'T';
end

Voltage = zeros(N_points,N_curves);
Current = zeros(N_points,N_curves);
Jph = zeros(1,N_curves);

AM1_5 = importdata(fullfile(Folder_name,'am15_300-1500.in'));

q = 1.60217662e-19;
AM1_5_wav = AM1_5.data(:,1);
AM1_5_flux = q*AM1_5.data(:,2);
for wav_i = 1:length(AM1_5_wav)-1
    AM1_5_flux(wav_i) = AM1_5_flux(wav_i)/(AM1_5_wav(wav_i+1)-AM1_5_wav(wav_i));
end

fig = figure;
if ~strcmp(Abs_name,'-'); fig.Position = fig.Position + [-200,0,400,0]; end
Lgd_text = cell(1,N_curves);
color_begin = [1,0,0]; %Red
color_end = [0,0,1]; %Blue
color = [linspace(color_begin(1),color_end(1),N_curves)',linspace(color_begin(2),color_end(2),N_curves)',linspace(color_begin(3),color_end(3),N_curves)'];
for i = 1:N_curves
    if strcmp(Type,'Irr'); T_i = 1; Irr_i = i; else; T_i = i; Irr_i = 1; end
    JV_Filename = fullfile(Folder_name,append(JV_name,'_',num2str(T_levels(T_i)),'deg_',num2str(Irr_levels(Irr_i)),'Wm2',Extension));
    Data = importdata(JV_Filename);
    Voltage(:,i) = Data.data(:,1);
    Current(:,i) = Data.data(:,2);
    if Flip_voltage; Voltage(:,i) = Voltage(:,i)*-1; end
    if Flip_current; Current(:,i) = Current(:,i)*-1; end

    if ~strcmp(Abs_name,'-'); subplot(1,2,1); end
    hold on; box on;
    plot(Voltage(:,i),Current(:,i),'color',color(i,:))

    if ~strcmp(Abs_name,'-')
        Abs_filename = fullfile(Folder_name,append(Abs_name,'_',num2str(T_levels(T_i)),'deg_',num2str(Irr_levels(Irr_i)),'Wm2',Extension));
        Data = importdata(Abs_filename);
        wav = Data.data(:,1);
        Abs = Data.data(:,Abs_ind);

        Jph(i) = trapz(wav,Abs.*AM1_5_flux)*Irr_levels(Irr_i)/1000;

        subplot(1,2,2)
        hold on; box on;
        plot(wav,Abs,'color',color(i,:))

    else
        Jph(i) = interp1(Voltage(:,i),Current(:,i),0);
    end
    if strcmp(Type,'Irr') 
        Lgd_text{i} = append('Irr = ',num2str(Irr_levels(i)),'W/m^2'); 
    else
        Lgd_text{i} = append('T = ',num2str(T_levels(i)),'C'); 
    end
end
if ~strcmp(Abs_name,'-'); subplot(1,2,1); end
xlabel('Voltage [V]')
ylabel('Current [A/m^2]')
xlim([0, max(max(Voltage))]);
ylim([0, 1.05*max(Current(1,:))]);
legend(Lgd_text,'location','northeast')

if ~strcmp(Abs_name,'-')
    subplot(1,2,2)
    xlabel('Wavelength [m]')
    ylabel('Absorption [-]')
    xlim([wav(1),wav(end)])
    ylim([0, 1]);
    legend(Lgd_text,'location','northeast')

end

drawnow

Ans_User =input('Are these curves correct? [Y/N]','s');

if strcmp(Ans_User,'Y')
    Check = 1;
else
    Check = 0;
end

end