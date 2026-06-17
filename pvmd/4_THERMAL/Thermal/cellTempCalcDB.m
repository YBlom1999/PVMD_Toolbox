function [T_cell] = cellTempCalcDB(T_amb, WS,Gm,settings)
%cellTempCalcDB Calculates the temperature of the cell with the
%Duffie-Beckman model.
% This function calculates the temperature of the cell with the
%Duffie-Beckman model.
%
% Parameters
% ----------
% T_amb : double
%   Ambient temperature
% WS : double
%   Windspeed
% Gm : double
%   incoming irradiance in the module
% settings : struct
%   the inputsettings of the thermal model
%
% Returns
% -------
% T_cell : double
%   The cell temperature
%
% Developed by A. Calcabrini
% Implemented in the Toolbox by Y. Blom
T_NOCT = settings.T_NOCT;
cell_eff = settings.cell_eff;

T_cell = T_amb+(T_NOCT-20)/800.*Gm*9.5./(5.7+3.8*WS)*(1-cell_eff/0.9);

end