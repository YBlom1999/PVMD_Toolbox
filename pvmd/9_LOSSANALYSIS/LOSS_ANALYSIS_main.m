function [LOSS_ANALYSIS_output, TOOLBOX_input] = LOSS_ANALYSIS_main(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_output)
%LOSS_ANALYSIS_MAIN Main file for the Loss analysis module in the PVMD toolbox
%
% This function calculates the losses that are present in the PV system
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Simulation results of the CELL module
% MODULE_output : struct
%   Simulation results of the MODULE module
% WEATHER_output : struct
%   Simulation results of the WEATHER module
% THERMAL_output : struct
%   Simulation results of the THERMAL module
% ELECTRIC_output : struct
%   Simulation results of the ELECTRIC module
% CONVERSION_output : struct
%   Simulation results of the CONVERSION module
%
% Returns
% -------
% LOSS_ANALYSIS_output : struct
%   Simulation results of the loss analysis module
% TOOLBOX_input : struct
%   Update of TOOLBOX_input by adding simulation parameters of loss analysis module
%
% Developed by Y. Blom

%---- Ask user input
if TOOLBOX_input.script == 1
    Run_Operating = TOOLBOX_input.LossAnalysis.Run_Operating;
    TrackingType = TOOLBOX_input.LossAnalysis.TrackingType;
    if ~isfield(TOOLBOX_input.LossAnalysis,'plotFigures')
        TOOLBOX_input.LossAnalysis.plotFigures = 0;
    end
else
    Answer = input("Do you want the loss analysis of operating conditions? [Y/N]","s");
    if Answer == 'Y'
        % Tracking should be improved before it can be used.
        %Answer = input("Do you want the simulate solar tracking? [Y/N]","s");
        %if Answer == 'Y'
        %    Q ={'Select type of tracking'};    %ask user
        %    TrackingType = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Azimuth tracking','Altitude tracking','Dual axis tracking'},'ListSize',[150,100]); %user choice
        %    Run_Operating = 0;
        %else
            Run_Operating = 1;
            TrackingType = 0;
        %end
    else
        Run_Operating = 0;
        TrackingType = 0;
    end

    TOOLBOX_input.LossAnalysis.Run_Operating = Run_Operating;
    TOOLBOX_input.LossAnalysis.Run_TrackingType = TrackingType;
    if strcmp(CELL_output.TYPE,'SHJ') || strcmp(CELL_output.TYPE,'BIF') || strcmp(CELL_output.TYPE,'T-F')
        A = get_user_input({'Bandgap energy of cell [eV]?'}, {'E_g'},{'1.12'});
        TOOLBOX_input.LossAnalysis.E_g = A.E_g;
    elseif strcmp(CELL_output.TYPE,'Tan') || strcmp(CELL_output.TYPE,'BIF-Tan')
        A = get_user_input({'Bandgap energy of top cell [eV]?','Bandgap energy of bottom cell [eV]?'}, {'E_g1','E_g2'},{'1.68','1.12'});
        TOOLBOX_input.LossAnalysis.E_g1 = A.E_g1;
        TOOLBOX_input.LossAnalysis.E_g2 = A.E_g2;
    end
    Answer = input("Do you want to plot the figures? [Y/N]","s");
    if Answer == 'Y'
        TOOLBOX_input.LossAnalysis.plotFigures = 1;
    else
        TOOLBOX_input.LossAnalysis.plotFigures = 0;
    end
end

%---- Calculate STC results
Losses_STC = Analysis_STC(TOOLBOX_input,CELL_output,MODULE_output,ELECTRIC_output,CONVERSION_output);
LOSS_ANALYSIS_output.Losses_STC = Losses_STC;


%---- Run loss analysis under operating conditions
if Run_Operating == 1
    disp('Operating conditions are simulated. This can take around 30 minutes.')
    Losses_Operating = Analysis_Operating(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_output);
    LOSS_ANALYSIS_output.Losses_Operating = Losses_Operating;
end
%---- Run loss analysis with tracking
if TrackingType ~= 0
    Losses_Tracking = Analysis_Tracking(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_out