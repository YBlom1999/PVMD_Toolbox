function [Tcell] = cellTempCalcSandia(T_amb,WS,Gm,thermal)
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

a = thermal.a;
b = thermal.b;

Tcell = Gm.*exp(a+b*WS)+T_amb;

end