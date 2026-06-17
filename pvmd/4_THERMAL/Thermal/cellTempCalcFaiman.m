function [T_cell] = cellTempCalcFaiman(T_amb,WS,Gm,settings)
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
% settings : struct
%   the input of the thermal model
%
% Returns
% -------
% T_cell : double
%   The cell temperature
%
% Developed by A. Calcabrini
% Implemented in the Toolbox by Y. Blom

U0 = settings.U0;
U1 = settings.U1;
alpha = settings.alpha;
cell_eff = settings.cell_eff;

T_cell = T_amb + Gm./(U0+U1*WS)*alpha*(1-cell_eff);


end