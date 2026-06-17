function [MODULE_output, WEATHER_output] = MODULE_main(TOOLBOX_input, CELL_output)
%WEATHER_MAIN Main file for the MODULE module in the PVMD toolbox
%
% This function calculates the sensitivity map of the PV system
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Output of the CELL block
%
% Returns
% -------
% MODULE_output : struct
%   Output of this module block
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by Karthik Ganapathi Subramanian, editted by Youri Blom

if TOOLBOX_input.runPeriodic && ~TOOLBOX_input.Scene.tracking %periodic simulation
    [MODULE_output] = MODULE_periodic(TOOLBOX_input, CELL_output);
    WEATHER_output = nan;
elseif TOOLBOX_input.runPeriodic && TOOLBOX_input.Scene.tracking %Tracking
    [WEATHER_output, MODULE_output] = Module_Tracking(TOOLBOX_input, CELL_output);
else %non-periodic simulation
    [MODULE_output] = MODULE_nonperiodic(TOOLBOX_input, CELL_output);
    WEATHER_output = nan;
end