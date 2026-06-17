function [specRef] = readSpectralReflectivity(materialName,wavelengths)
% readSpectralReflectivity reads the spectral reflectivity of a material
%
% This function uses the specified filename and interpolates it for the
% specified wavelenghts
%
% Parameters
% ----------
% materialName : string
%   The material for which the reflection is obtained
% wavelenghts: double
%   The wavelengths for which the spectral reflectivity should be
%   calculated
%
% Returns
% -------
% specRef : double
%   The spectral reflectivity of the material at the specified wavelengths
%
% Developed by Youri Blom

T = readtable([materialName,'.xlsx']);
wavOriginal = T{:,1};
reflectionOriginal = T{:,2};

specRef = interp1(wavOriginal,reflectionOriginal,wavelengths*1e9);

end