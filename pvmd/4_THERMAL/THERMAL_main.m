function [THERMAL_output, TOOLBOX_input] = THERMAL_main(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output)
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
% TOOLBOX_input : struct
%   Update of TOOLBOX_input by adding simulation parameters of thermal module
% 
% Developed by unknown (A. Jamodkar? E. Garcia?). Improved by A. Nour.
% Commented by A. Alcaniz
% Refactored by M. Kok



%---- Set user parameters if missing
TOOLBOX_input = set_config_thermal(TOOLBOX_input);

    

%---- Calculate cell temperatures

if TOOLBOX_input.thermal.Type == 1
    THERMAL_output = calculate_cell_temp(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output);
elseif TOOLBOX_input.thermal.Type == 2
    [THERMAL_output, TOOLBOX_input] = PVT_col_main(TOOLBOX_input, MODULE_output, WEATHER_output);
end


%---- Save THERMAL output to disk 
if ~isempty(THERMAL_output) && exist('TOOLBOX_input.save','var') && ~isempty(TOOLBOX_input.save.folder)
    disp('Saving thermal output...')
    file = fullfile(TOOLBOX_input.save.folder, 'THERMAL_output.mat');
    save(file, 'THERMAL_output')
end
end
