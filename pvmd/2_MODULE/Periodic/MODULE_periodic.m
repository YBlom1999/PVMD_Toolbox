function [MODULE_output] = MODULE_periodic(TOOLBOX_input, CELL_output)
%WEATHER_periodic Calculates the sensitivity map for periodic simulations
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
% Developed by Karthik Ganapathi Subramanian


%load input values
TOOLBOX_input.runPeriodic = true; %setting periodic simulation as true
TOOLBOX_input.Scene.NonPeriodic_Environment = '';

if TOOLBOX_input.Scene.loadSimulation
    load(TOOLBOX_input.Scene.LuxLoadFile,'MODULE_output')
    
elseif TOOLBOX_input.Scene.runLUX
    [MODULE_output] = runLUX_main(TOOLBOX_input, CELL_output);
elseif TOOLBOX_input.Scene.runBackwardTracer
    [MODULE_output] = runBackwardTracer_main(TOOLBOX_input, CELL_output);
end

end
