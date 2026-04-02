function [pvt_collector_output, TOOLBOX_input] = PVT_col_main(TOOLBOX_input, MODULE_output, WEATHER_output)
% PVT_col_MAIN Main file for the PVT module in the PVMD toolbox
%
%This function calculates the thermal & PV yield of PVT collector
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% MODULE_output : struct
%   Output of this module block
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% pvt_collector_output : struct
%   Simulation results of the PVT module
% TOOLBOX_input : struct
%   Update of TOOLBOX_input by adding simulation parameters of thermal module
% 
% Developed by ZUA.

% Load the chosen function
fh = str2func(TOOLBOX_input.thermal.runpvt(1:end-2)); % Get the function handle

% Call the chosen function and get the output
[pvt_collector_output] = fh(MODULE_output, WEATHER_output);

%% Plots
% Plots the energy yield and irradiance
TOOLBOX_input.runPeriodic;
plot_EnergyYield_PVT(pvt_collector_output,WEATHER_output);

% Adjust the temperature array for electrical simulation % Y.
pvt_collector_output.T = repmat(pvt_collector_output.T',[1,72]);

if ~TOOLBOX_input.script
    if pvt_collector_output.total_electrical_output>0
        disp('PVT Collector calculation finished, PVT yield calculated.');
        disp(['Total Thermal Yield in given period:', num2str(pvt_collector_output.total_thermal_output), 'kWh/m^2']);
        disp(['Total Electrical Yield in given period:', num2str(pvt_collector_output.total_electrical_output), 'kWh/m^2']);
    elseif pvt_collector_output.total_thermal_output==0
        disp('Heat Pump Calculation finnished.');
    else
        disp('Thermal Collector calculation finished, Thermal yield calculated.');
        disp(['Total Thermal Yield in given period:', num2str(pvt_collector_output.total_thermal_output), 'kWh/m^2']);
    end
end
end
