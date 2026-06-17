function [LOSS_ANALYSIS_output] = LOSS_ANALYSIS_main(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_output)
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


Run_Operating = TOOLBOX_input.LossAnalysis.Run_Operating;
TrackingType = TOOLBOX_input.LossAnalysis.TrackingType;
if ~isfield(TOOLBOX_input.LossAnalysis,'plotFigures')
    TOOLBOX_input.LossAnalysis.plotFigures = 0;
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
    Losses_Tracking = Analysis_Tracking(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_output,TrackingType);
    LOSS_ANALYSIS_output.Losses_Tracking = Losses_Tracking;
end

end