function [Tcell] = cellTempCalcFaiman(T_amb,WS,Gm,thermal)
%cellTempCalcFaiman Calculates the temperature of the cell with the
%Faiman model.
% This function calculates the temperature of the cell with the
%Faiman model.
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

U0 = thermal.U0;
U1 = thermal.U1;
alpha = thermal.alpha;
cell_eff = thermal.cell_eff;

Tcell = T_amb + Gm./(U0+U1*WS)*alpha*(1-cell_eff);


end