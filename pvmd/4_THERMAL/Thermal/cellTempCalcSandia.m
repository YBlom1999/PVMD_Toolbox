function [T_cell] = cellTempCalcSandia(T_amb,WS,Gm,settings)
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

a = settings.a;
b = settings.b;

T_cell = Gm.*exp(a+b*WS)+T_amb;

end