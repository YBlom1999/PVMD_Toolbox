function [MODULE_output, TOOLBOX_input, WEATHER_output] = MODULE_main(TOOLBOX_input, CELL_output)
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

if TOOLBOX_input.script == 0
    usrops = {'Module with Periodic Properties','Two-axis tracking','Non-periodic Simulation'};
    promptstring = 'Select Operation';
    idx = listdlg('PromptString',promptstring...
        ,'ListString',usrops,'SelectionMode','single','ListSize',[200 100]);
    if idx == 1
        TOOLBOX_input.runPeriodic = true;
        TOOLBOX_input.Scene.tracking = false;
    elseif idx == 2
        TOOLBOX_input.runPeriodic = true;
        TOOLBOX_input.Scene.tracking = true;
    else
        TOOLBOX_input.runPeriodic = false;
        TOOLBOX_input.Scene.tracking = false;
    end
    if isempty(idx)
        cd ..
        return
    end
end

if TOOLBOX_input.runPeriodic && ~TOOLBOX_input.Scene.tracking %periodic simulation
    [MODULE_output, TOOLBOX_input] = MODULE_periodic(TOOLBOX_input, CELL_output);
    WEATHER_output = nan;
elseif TOOLBOX_input.runPeriodic && TOOLBOX_input.Scene.tracking %Two axis tracking
    [WEATHER_output, MODULE_output, TOOLBOX_input] = TwoAxisTracking(TOOLBOX_input, CELL_output);
else %non-periodic simulation
    [MODULE_output, TOOLBOX_input] = MODULE_nonperiodic(TOOLBOX_input, CELL_output);
    WEATHER_output = nan;
end