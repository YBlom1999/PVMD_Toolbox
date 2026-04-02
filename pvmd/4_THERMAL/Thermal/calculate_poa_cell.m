function poa_cell = calculate_poa_cell(periodic, WEATHER)
%CALCULATE_POA_CELL Calculate plane of array irradiance in the cell
%
% Parameters
% ----------
% periodic : logical
%   Indicates if the simulation is periodic or not
% num_modules : double
%   Number of modules in the system. 1 for periodic simulations
% WEATHER : struct
%   Simulation results of the weather block
%
% Returns
% -------
% poa_cell : cell
%   Plane of array irradiance per cell and per module

if ~periodic
    poa_cell = WEATHER.A;      
else
    poa_cell = {WEATHER.A};
end

end
