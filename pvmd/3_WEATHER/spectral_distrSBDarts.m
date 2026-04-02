function [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,air_mass,wav_new] = spectral_distrSBDarts(wav)
%spectral_distrSBDarts Calculate the relative spectral distribution with
%SBDarts
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
% RSD_i_dir : double
%   Relative spectral distribution in terms of irradiance for direct irradiation [W/m2]
% RSD_f_dir : double
%   Relative spectral distribution in terms of photon flux for direct irradiation [#/m2]
% RSD_i_dif : double
%   Relative spectral distribution in terms of irradiance for diffuse irradiation [W/m2]
% RSD_f_dif : double
%   Relative spectral distribution in terms of photon flux for diffuse irradiation [#/m2]
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
load('constants/spectraSBDART.mat','sbdartSpec')
air_mass = 1./sind(90-sbdartSpec.solarZenith(1:end-1));
wavelengths = sbdartSpec.lambda/1e3; %Wavelength in um
initial_wav = wavelengths(1);
final_wav = wavelengths(end);
spectr_power_dens_dir = sbdartSpec.direct(:,1:end-1,:);
spectr_power_dens_dif = sbdartSpec.diffuse(:,1:end-1,:);

% Prepare wavelength boundaries, adding the extreme wavelengths (280 and
% 4000 nm) if needed
stepsize = mean(wav(2:end)-wav(1:end-1));
wav_new = wav;
if wav_new(1)> initial_wav, wav_new = [initial_wav:stepsize:wav_new(1),wav_new(2:end)']'; end
if wav_new(end)< final_wav, wav_new = [wav_new(1:end-1)',wav_new(end):stepsize:final_wav]'; end

% Interpolate the spectrum at the new wavelength range
RSD_i_dir = interp1(wavelengths,spectr_power_dens_dir,wav_new);
RSD_i_dif = interp1(wavelengths,spectr_power_dens_dif,wav_new);

% Spectral photon flux density [photons/m2/s/nm]
RSD_f_dir = RSD_i_dir.*wav_new*1e-6/(h*c);
RSD_f_dif = RSD_i_dif.*wav_new*1e-6/(h*c);


% Perform normalization by dividing the spectral distributions over the
% total irradiance
total_irr_dir = trapz(wav_new,RSD_i_dir);
total_irr_dir(total_irr_dir == 0) = 1;
RSD_i_dir = RSD_i_dir./total_irr_dir;
RSD_f_dir = RSD_f_dir./total_irr_dir;

total_irr_dif = trapz(wav_new,RSD_i_dif);
RSD_i_dif = RSD_i_dif./total_irr_dif;
RSD_f_dif = RSD_f_dif./total_irr_dif;
end
