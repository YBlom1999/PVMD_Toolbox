function [pvModule,weather] = prepare_incropera(TOOLBOX_input, WEATHER, CELL_output,cell_i)
% prepare_incropera prepares the correct inputs for the incropera model.
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% WEATHER : struct
%   Results from the WEATHER module
% CELL_output : struct
%   Results from the CELL module
% cell_i: double
%   Index of which cell is considered
%
% Returns
% -------
% pvModule : struct
%   Properties of the simulated module
% weather: struct
%   The relevant weather parameters
%
% Developed by A. Calcabrini, commented by Y. Blom

pvModule.tilt = TOOLBOX_input.Scene.module_mounting.ModTilt;
pvModule.azimuth = TOOLBOX_input.Scene.module_mounting.ModAzimuth;
pvModule.albedo = TOOLBOX_input.Scene.module_mounting.Albedo;

CR = TOOLBOX_input.Scene.module_mounting.CellRows;
CC = TOOLBOX_input.Scene.module_mounting.CellColumns;
CS = TOOLBOX_input.Scene.module_mounting.CellSpacing*1e-2;
ES = TOOLBOX_input.Scene.module_mounting.EdgeSpacing*1e-2;
CL = TOOLBOX_input.Scene.module_mounting.CellLength*1e-2;
CW = TOOLBOX_input.Scene.module_mounting.CellWidth*1e-2;
pvModule.length = CR*CL + (CR-1)*CS + 2*ES;
pvModule.width = CC*CW + (CC-1)*CS + 2*ES;
pvModule.layers = TOOLBOX_input.thermal.layers;


[weather_data, ~] = load_meteonorm_data(TOOLBOX_input);
weather.dhi = weather_data(:,8);
weather.dni = weather_data(:,7);
weather.ghi = weather_data(:,11);
weather.tAmb = weather_data(:,9);
weather.windSpeed = weather_data(:,10);
weather.solarAzi = weather_data(:,5);
weather.solarAlti = weather_data(:,6);
weather.tGnd = (weather.tAmb+273.15)+(0.015*(1-pvModule.albedo).*weather.ghi-0.7).*exp(-0.09*weather.windSpeed)-273.15;%(Celsius)
weather.tSky = 0.0552*(weather.tAmb+273.15).^1.5-273.15;%(Celsius) Obtained from the Solar Energy Book
qGen = zeros(TOOLBOX_input.thermal.Nlayers,length(WEATHER.Irr));
assignment = TOOLBOX_input.thermal.assignment_layers;
for i = 2:length(assignment)-1
    if sum(i == CELL_output.CELL_FRONT.Absmat_ind)
        qGen(assignment(i),:) = qGen(assignment(i),:)+WEATHER.A(:,cell_i,i-1)'*(1-TOOLBOX_input.thermal.Efficiency/100);
    else
        qGen(assignment(i),:) = qGen(assignment(i),:)+WEATHER.A(:,cell_i,i-1)';
    end
end
weather.qGen = qGen;
end