function [THERMAL_output] = THERMAL_main(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output)
%THERMAL_MAIN Main file for the Thermal module in the PVMD toolbox
%
% This function calculates the temperature of the cells in the system
% considering the weather and the absorbed irradiance
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% THERMAL_output : struct
%   Simulation results of the thermal module
% 
% Developed by unknown (A. Jamodkar? E. Garcia?). Improved by A. Nour.
% Commented by A. Alcaniz
% Refactored by M. Kok

%---- Calculate cell temperatures
if TOOLBOX_input.thermal.Type == 1
    THERMAL_output = calculate_cell_temp(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output);
elseif TOOLBOX_input.thermal.Type == 2
    THERMAL_output = PVT_col_main(TOOLBOX_input, MODULE_output, WEATHER_output);
end
end