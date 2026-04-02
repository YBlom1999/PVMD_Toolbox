function rec_irradiance= surrogate_fun_per_hour(x, CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
    albedo,Bi)
%surrogate_fun_per_hour Is the objective function for the surrogate
%optimization
%
% This function calculates received irradiance for a given orientation
%
% Parameters
% ----------
% x : double
%   Parameters that are optimized, x consists of [azim, tilt]
% CELL_output : struct
%   Output of the CELL block
% azim_vertex : double
%   the azimuth of all vertices
% zenith_vertex : double
%   the zenith of all vertices
% wav_WEATHER : double
%   the wavelengths that are considered by the WEATHER module
% albedo : double
%   the spectral albedo of the ground
% Bi : double
%   the irradiance from all vertices
%
% Returns
% -------
% rec_irradiance : double
%   The received irradiance for this orientation
%
% Developed by Orestis Chatzilampos, integrated by Youri Blom
wav=wav_WEATHER; 
azim=x(1);
tilt=x(2);


SM=calculate_SM_new_single_orientation(CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
    albedo,azim,tilt);

rec_irradiance=trapz(wav,sum(Bi.*squeeze(SM(:,end,:)))');

end