function [Tcell] = cellTempCalcDB(T_amb, WS,Gm,thermal)
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
% thermal : struct
%   the input of the thermal model
%
% Returns
% -------
% Tcell : double
%   The cell temperature
%
% Developed by A. Calcabrini
% Implemented in the Toolbox by Y. Blom
T_NOCT = thermal.T_NOCT;
cell_eff = thermal.cell_eff;

Tcell = T_amb+(T_NOCT-20)/800.*Gm*9.5./(5.7+3.8*WS)*(1-cell_eff/0.9);

end