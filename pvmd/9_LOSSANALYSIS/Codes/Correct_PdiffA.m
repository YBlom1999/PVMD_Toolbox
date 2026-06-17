function [P_diffA, I_abs] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs,T_cell,V_opt,IVtype)
% Correct_PdiffA updates the parasitc absorption losses and the absorbed
% current
%
% This function accounts for the fact that the Jsc does not exactly equal
% the absorbed current, and corrects for this in the parasitic absorption
% losses and absorbed current
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% MODULE_output : struct
%   Simulation results of the MODULE module
% P_diffA : double
%   The original losses due to the Parasitic absorption
% I_abs : double
%   The original absorbed current
% T_cell : double
%   The temperature of the cell
% Vopt : double
%   The optimal voltage under ideal conditions
% IVtype : cell
%   The filename of the IV characteristics
%
% Returns
% -------
% P_diffA : double
%   The updated losses due to the Parasitic absorption
% I_abs : double
%   The updated absorbed current
%
% Developed by Y. Blom

%Properties of the module
N_cells = MODULE_output.N;
A_cell = MODULE_output.A;
Sf = TOOLBOX_input.electric.shading/100;
A_eff = N_cells*A_cell*(1-Sf);
%If the Isc is larger than the absorbed current (possible for SHJ), this
%needs to be taken from the parasitic absorption.

[paramForwardBias,~,~] = readIVcurveFile(IVtype);
Jph_J = paramForwardBias.Jph_J;
Jph_T = paramForwardBias.Jph_T;
I_sc=(A_eff/N_cells)*polyval(Jph_T,T_cell).*polyval(Jph_J,I_abs/(A_eff/N_cells))./polyval(Jph_T,298.15);
Change_ParAbs = I_abs < I_sc;
P_diffA(Change_ParAbs) = P_diffA(Change_ParAbs)-N_cells*V_opt(Change_ParAbs).*(I_sc(Change_ParAbs)-I_abs(Change_ParAbs));
I_abs(Change_ParAbs) = I_sc(Change_ParAbs);
end