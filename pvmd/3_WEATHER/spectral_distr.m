function [RSD_i,RSD_f,air_mass] = spectral_distr(wav)
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
%
% Developed by unknown (V. Muthukumar? E. Garcia?). Commented by A.
% Alcaniz


% Load constants and air mass spectra (spectral irradiance for 2002
% wavelengths, AM1.0, AM1.5, etc)
load('constants/weather_params.mat', 'h','c')
load('constants/spectral_data.mat','air_mass','wavelengths',...
    'initial_wav','final_wav','spectr_power_dens')

% Prepare wavelength boundaries, adding the extreme wavelengths (280 and
% 4000 nm) if needed
wav_b = (wav(1:end-1)+wav(2:end))/2;
wav_b = [wav(1)-(wav(2)-wav(1))/2;wav_b;wav(end)+(wav(end)-wav(end-1))/2];
if wav_b(1)> initial_wav, wav_b = [initial_wav;wav_b]; end
if wav_b(end)< final_wav, wav_b = [wav_b;final_wav]; end

% Spectral photon flux density [photons/m2/s/nm]
spectr_photon_flux = spectr_power_dens.*wavelengths*1e-6/(h*c);

% For every wavelength integrate the irradiance and the spectral photon
% flux
number_wavelengths = length(wav_b)-1;  
RSD_i = zeros(number_wavelengths,length(air_mass));
RSD_f = zeros(size(RSD_i));
for w = 1:number_wavelengths
    ix = wavelengths>wav_b(w)-0.0001 & wavelengths<=wav_b(w+1)+0.0001;
    RSD_i(w,:) = trapz(wavelengths(ix),spectr_power_dens(ix,:));
    RSD_f(w,:) = trapz(wavelengths(ix),spectr_photon_flux(ix,:));
end

% Perform normalization by dividing the spectral distributions over the
% total irradiance
total_irr = ones(number_wavelengths,1)*trapz(wavelengths,spectr_power_dens);
RSD_i = RSD_i./total_irr;
RSD_f = RSD_f./total_irr;

end
