function [wl_um,j,status] = read_spectrum(SpecName,root)
% Reads spectrum from sheet in spectrum.xlsx file and converts to current.

% ---INPUT---
% name      spectrum name e.g. 'AM1_5g' (should match sheet name in .xlsx)

% ---OUTPUT---
% wl_um     interval midpoint wavlength [um]
% j         interval implied photocurrent [mA/cm2]
% status    0 = ok, 1 = data not found, 2 = data format incorrect

%---constants of nature---
h =  6.62607015e-34;    % Planck constant [J s]
c =  299792458;         % speed of light [m/s]
q =  1.60217663e-19;    % elementary charge [C]
%---

SpecFilePath = fullfile(root,'spectrum.xlsx');
sheets = sheetnames(SpecFilePath);   % check sheets in the xlsx file

if any(sheets == SpecName)                  % if sheet name matches
    spectrum = readtable(SpecFilePath,"Sheet",SpecName); % read data

    wl_unit = spectrum.Properties.VariableNames{1};   % wavelength
    spec_unit = spectrum.Properties.VariableNames{2}; % intensity

    % wavelength must be in nm, intensity must be in W/m2/nm
    % TODO: could be made to accept um,m,eV (just like read_nk)
    if strcmp(wl_unit,'nm') && strcmp(spec_unit,'W_m2_nm')

        % wavelength interval midpoints [um]
        wl_um = (spectrum.nm(1:end-1)+spectrum.nm(2:end))/2000;
        E_ph = 1e6*h*c./wl_um;       % corresponding photon energy [J]

        % spectral intensity [W/m2/nm] at those midpoints
        intensity = (spectrum.W_m2_nm(1:end-1)+spectrum.W_m2_nm(2:end))/2;

        % wavelength interval widths [nm]
        width = spectrum.nm(2:end) - spectrum.nm(1:end-1);

        % nr of photons in interval [photons/m2/s]
        nr_photons = intensity .* width ./ E_ph;

        % implied photocurrent of photons in the interval [mA/cm2]
        j = nr_photons * q /10;     %x1000 A-->mA, /10000 m2-->cm2

        status = 0;
    else                                % if unit are not correct 
        wl_um = [];
        j = [];
        status = 2;
    end


else                                    % if no matching sheet found
    wl_um = [];
    j = [];
    status = 1;
end
