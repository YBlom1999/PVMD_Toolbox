function [Lay, Int, absmat, rback, Type, Submod_ind, SET, wav] = prepareInputCell(TOOLBOX_input)
% prepareInputCell prepares the input of the cell simulation
%
% This function reads the GenPro file and extracts all relevant information
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
%
% Returns
% -------
% Lay: struct
%   The layers of the simulated cell
% Int: struct
%   The interfaces between the layers
% Absmat: double
%   Index of which layers or interfaces are absorber materials
% rback: double
%   Reflectance of backsheet (nan for bifacial modules)
% Type: string
%   The technology type that is used for the module
% Submod_ind: double
%   Index of which absorber materials belong to which submodule
% SET: struct
%   An overview of the settings that are used
% wav: double
%   The wavelengths that are simulated
%
% Developed by Youri Blom

%---set angles of incidence---
SET = read_settingsGP(TOOLBOX_input.settings);                         %load GenPro4 settings
SET.save_ang_tree = 0;
SET.results_to_file = 0;
SET.plot_fig_input = 0;
SET.plot_fig_output = 0;

%get info from seperately defined GenProFile
[DeviceTable, DeviceInfo] = read_device(TOOLBOX_input.deviceOptic.GenProFile);
[Lay,Int,~,absmat] = table2struct_GP(DeviceTable,SET);

rback = str2double(DeviceInfo{1,2});
Type = DeviceInfo{2,2};
Submod_ind = str2num(DeviceInfo{3,2});

wav = (SET.start_wavelength_nm:SET.step_wavelength_nm:SET.stop_wavelength_nm)/1e9;

if isnan(rback)                   %if rback is empty, bificial module is used
    SET.skip_dark = 0;                               %disables an acceleration trick that doesn't work for rear side illumination
end

end