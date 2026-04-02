function [RSD_i,RSD_f,air_mass,wav_new] = spectral_distrYB(wav)
%SPECTRAL_DISTR Calculate the relative spectral distribution
%
% Convert all spectra to relative spectral distibution over generated
% wavelength boundaries. Two types:
% 1) irradiance in interval, divided by total irradiance [-]
% 2) photon flux in interval, divided by total irradiance [photons/J]
%
% Parameters
% ----------
% wav : double
%   Wavelengths [um] interval centers (might not cover entire spectrum)
% 
% Returns
% -------
% RSDi : double
%   Relative spectral distribution in terms of irradiance [W/m2]
% RSDf : double
%   Relative spectral distribution in terms of photon flux [W/m2]
% AM : double
%   Air mass numbers for which RSD is calculated
% wav_new: double
%   Wavelength range for which the RSDi and RSDf are calculated [um]
%
% Developed by unknown (V. Muthukumar? E. Garcia?). %Changed by Y, Blom. %Commented by A.
% Alcaniz


% Load constants and air mass spectra (spectral irradiance for 2002
% wavelengths, AM1.0, AM1.5, etc)
load('constants/weather_params.mat', 'h','c')
load('constants/spectral_data.mat','air_mass','wavelengths',...
    'initial_wav','final_wav','spectr_power_dens')

% Prepare wavelength boundaries, adding the extreme wavelengths (280 and
% 4000 nm) if needed
stepsize = mean(wav(2:end)-wav(1:end-1));
wav_new = wav;
if wav_new(1)> initial_wav, wav_new = [initial_wav:stepsize:wav_new(1),wav_new(2:end)']'; end
if wav_new(end)< final_wav, wav_new = [wav_new(1:end-1)',wav_new(end):stepsize:final_wav]'; end

% Interpolate the spectrum at the new wavelength range
RSD_i = interp1(wavelengths,spectr_power_dens,wav_new);

% Spectral photon flux density [photons/m2/s/nm]
RSD_f = RSD_i.*wav_new*1e-6/(h*c);


% Perform normalization by dividing the spectral distributions over the
% total irradiance
total_irr = trapz(wav_new,RSD_i);
RSD_i = RSD_i./total_irr;
RSD_f = RSD_f./total_irr;


end
