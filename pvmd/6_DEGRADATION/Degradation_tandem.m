function [CELL_output_new,ELECTRIC_output_new,k_dis,k_mois,k_LID,k_TC,Time,Rcon_new] = Degradation_tandem(TOOLBOX_input,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONSTANTS)
%Degradation_single file for the degradation module of tandem modules
%
% This function calculates the degradation rate of a tandem module
% 
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% MODULE_output : struct
%   Simulation results of the MODULE module
% WEATHER_output : struct
%   Simulation results of the WEATHER module
% THERMAL_output : struct
%   Simulation results of the THERMAL module
% ELECTRIC_output : struct
%   Simulation results of the ELECTRIC module
% Parameters : double
%   The parameters of the degradation mechanisms
% CONSTANTS : double
%   Physical constants
%
% Returns
% -------
% CELL_output_new : struct
%   New simulation results of the CELL module
% ELECTRIC_output_new : struct
%   New simulation results of the ELECTRIC module
% k_dis: double
%   The degradation rate due to discoloration
% k_mois: double
%   The degradation rate due to moisture ingress
% k_LID: double
%   The degradation rate due to light induced degradation
% k_TC: double
%   The degradation rate due to thermal cycling
% Time: double
%   The time for which degradation is calculated
% Rcon_new:
%   The new value for the interconnection resistance
%
% Developed by by Youri Blom

%Degradation parameters
A_mois = TOOLBOX_input.Degradation.A_mois;
C_mois = TOOLBOX_input.Degradation.C_mois;
Ea_mois = TOOLBOX_input.Degradation.Ea_mois;
A_dis = TOOLBOX_input.Degradation.A_dis;
Ea_dis = TOOLBOX_input.Degradation.Ea_dis;
factor_max = TOOLBOX_input.Degradation.factor_max;
C_LID = TOOLBOX_input.Degradation.C_LID;
A_TC = TOOLBOX_input.Degradation.A_TC;
Ea_TC = TOOLBOX_input.Degradation.Ea_TC;
n_TC = TOOLBOX_input.Degradation.n_TC;
b_TC = TOOLBOX_input.Degradation.b_TC;

%Module parameters
numCells = MODULE_output.N;
Rcon = TOOLBOX_input.electric.resistance;

N_years = TOOLBOX_input.Degradation.N_years;

%Make repeating arrays of the stress factors
T_repeat = repmat(mean(THERMAL_output.T{1},2)+273.15,[N_years,1])';
UV_repeat = repmat(mean(WEATHER_output.UV{1},2),[N_years,1])';
Jph_repeat = repmat(CONSTANTS.q*mean(WEATHER_output.J{1}(:,:,2),2),[N_years,1])';

%Calculate degradation rates
Time = 1:N_years*8760;
[k_dis,CELL_output_new,MODULE_output_new,WEATHER_output_new,THERMAL_output_new] = DegDiscoloration(TOOLBOX_input,T_repeat,UV_repeat,A_dis,Ea_dis);
[k_mois] = DegMoistureIngress(TOOLBOX_input,T_repeat,Ea_mois,A_mois,C_mois,Time);
[k_LID,factor_LID] = DegLID(ELECTRIC_output,Jph_repeat,factor_max,C_LID,Time,TOOLBOX_input.constants);
[k_TC,Rcon_new] = DegThermalCycling(ELECTRIC_output,T_repeat,A_TC,Ea_TC,n_TC,b_TC,Rcon,numCells,TOOLBOX_input.constants);

TOOLBOX_input_new = TOOLBOX_input;

%Setting for bottom cell
TOOLBOX_input_new.electric.deg_sil.k_mois = sum(k_mois);
TOOLBOX_input_new.electric.deg_sil.factor_LID = factor_LID;
TOOLBOX_input_new.electric.deg_sil.filename = TOOLBOX_input.Degradation.MoistureDeg_file;
TOOLBOX_input_new.electric.resistance = Rcon_new;


%Setting for top cell
TOOLBOX_input_new.electric.k_needed{1} = N_years*TOOLBOX_input.Degradation.FactorDegPerov/100;
TOOLBOX_input_new.electric.degScenario{1} = TOOLBOX_input.Degradation.DegScenario;

%Run new electrical simulation
TOOLBOX_input_new.script = 1;
TOOLBOX_input_new.electric.electricplot = 0;
[ELECTRIC_output_new] = ELECTRIC_main(TOOLBOX_input_new, CELL_output_new, MODULE_output_new,WEATHER_output_new,THERMAL_output_new);

end