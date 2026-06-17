function [DEGRADATION_output] = DEGRADATION_main(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output)
% DEGRADATION_MAIN Main file for the degradation module in the PVMD toolbox
%
% This function calculates the degradation rate of a PV module under
% outdoors operating conditions
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
%
% Returns
% -------
% DEGRADATION_output : struct
%   Simulation results of the degradation module
%
% Developed by by Youri Blom


TOOLBOX_input.Degradation = readDegradationSettings(TOOLBOX_input.Degradation,TOOLBOX_input.settings);

if strcmp(CELL_output.TYPE,'Tan')
    [CELL_output_new,ELECTRIC_output_new,k_dis,k_mois,k_LID,k_TC,Time,Rcon_new] = Degradation_tandem(TOOLBOX_input,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,TOOLBOX_input.constants);
else
    [CELL_output_new,ELECTRIC_output_new,k_dis,k_mois,k_LID,k_TC,Time,Rcon_new] = Degradation_single(TOOLBOX_input,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,TOOLBOX_input.constants);
end




%Total Degradation
k_total = k_dis + k_mois+k_LID+k_TC;

numCells = MODULE_output.N;
Rcon = TOOLBOX_input.electric.resistance;
% plotDegradation(CELL_output,CELL_output_new,ELECTRIC_output,ELECTRIC_output_new,Time,k_dis,k_mois,k_LID,k_TC,k_total,numCells,Rcon,Rcon_new,CONSTANTS);

DEGRADATION_output.Time = Time;
DEGRADATION_output.k_mois = k_mois;
DEGRADATION_output.k_dis = k_dis;
DEGRADATION_output.k_LID = k_LID;
DEGRADATION_output.k_TC = k_TC;
DEGRADATION_output.k_total = k_total;
DEGRADATION_output.CELL_output_new = CELL_output_new;
DEGRADATION_output.ELECTRIC_output_new = ELECTRIC_output_new;
DEGRADATION_output.Rcon_new = Rcon_new;

end