function [DEGRADATION_output, TOOLBOX_input] = DEGRADATION_main(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output)
%DEGRADATION_MAIN Main file for the degradation module in the PVMD toolbox
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
% TOOLBOX_input : struct
%   Update of TOOLBOX_input by adding simulation parameters of degradation module
%
% Developed by by Youri Blom

%constants
CONSTANTS.k_b = 1.380649e-23;
CONSTANTS.q = 1.60217662e-19;
CONSTANTS.T_STC = 298.15;

Type = CELL_output.TYPE;

if ~TOOLBOX_input.script
    [TOOLBOX_input] = getDegradationInputs(TOOLBOX_input);
end

Parameters = readDegradationParameters(TOOLBOX_input);

if strcmp(Type,'Tan')
    [CELL_output_new,ELECTRIC_output_new,k_dis,k_mois,k_LID,k_TC,Time,Rcon_new] = Degradation_tandem(TOOLBOX_input,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,Parameters,CONSTANTS);
else
    [CELL_output_new,ELECTRIC_output_new,k_dis,k_mois,k_LID,k_TC,Time,Rcon_new] = Degradation_single(TOOLBOX_input,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,Parameters,CONSTANTS);
end




%Total Degradation
k_total = k_dis + k_mois+k_LID+k_TC;

numCells = MODULE_output.N;
Rcon = TOOLBOX_input.electric.resistance;
% plotDegradation(CELL_output,CELL_output_new,ELECTRIC_output,ELECTRIC_output_new,Time,k_dis,k_mois,k_LID,k_TC,k_total,numCells,Rcon,Rcon_new,CONSTANTS);


DEGRADATION_output.k_mois = k_mois;
DEGRADATION_output.k_dis = k_dis;
DEGRADATION_output.k_LID = k_LID;
DEGRADATION_output.k_TC = k_TC;
DEGRADATION_output.k_total = k_total;
DEGRADATION_output.CELL_output_new = CELL_output_new;
DEGRADATION_output.ELECTRIC_output_new = ELECTRIC_output_new;
DEGRADATION_output.Rcon_new = Rcon_new;

end