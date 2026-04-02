function output_irradiance = modify_irradiance(TOOLBOX_input, input_irradiance)
%MODIFY_IRRADIANCE Modify the irradiance with the PV system surroundings
%
% LiDAR data is employed to reconstruct the horizon at the PV module
% position. Then the SVF and the shading casted by objects are calculated
% and used to modify the measured irradiance at the weather station.
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% input_irradiance : double
%   Input weather data for selected time period. Columns in order: sun
%   azimuth, sun altitude, DNI, DHI
% var_argin : ??
%   Optional arguments used when running inhomogenous simulations
%
% Returns
% -------
% DNI_mod : double
%   DNI corrected by the shading factor
% DHI_mod : double
%   DHI corrected by the SVF
%
% Developed by A. Alcaniz

% Calculate the skyline and SVF if applicable
if TOOLBOX_input.irradiation.recalculate_horizon
    disp('Reconstructing the horizon...')
    horizon_reconstruction(TOOLBOX_input)
    disp('Horizon reconstructed')
end

% Load the reconstructed horizon
[~,~,data_folder] = get_folder_structure;
rec_hor_dir = fullfile(data_folder,'Weather','Reconstructed Horizons');
skyline_file = fullfile(rec_hor_dir,TOOLBOX_input.irradiation.skyline_file);
load(skyline_file,'skyline','SVF')

% Plot the skyline if indicated
if TOOLBOX_input.irradiation.plot_skyline
    plot_skyline(skyline)
end

% Compute the shading factor
shading_factor = compute_shading_factor(...
    skyline, ...
    input_irradiance(:,1) + 180, ...
    input_irradiance(:,2));

% Create the variable where the results will be stored
output_irradiance = input_irradiance(:,3:4);

% Reduce DNI by the shading factor
output_irradiance(:,1) = input_irradiance(:,3).*shading_factor;

% Reduce DHI by the sky view factor (SVF)
output_irradiance(:,2) = input_irradiance(:,4)*SVF;
end


function plot_skyline(skyline)
%PLOT_SKYLINE plot the skyline
%
% Parameters
% ----------
% skyline : cell/double
%   Skyline of the selected location

figure(); box on; hold off
set(gcf, 'Position',  [200 300 1200 400])
if iscell(skyline)
    if all(size(skyline) == [1,2])
        x_azim = skyline{1,1}; y_alti = skyline{1,2};
    elseif all(size(skyline) == [1,1])
        skl = skyline{1,1};
        x_azim = skl(:,1); y_alti = skl(:,2);
    end
    area(x_azim,y_alti,'FaceColor',[0.7 0.7 0.7])
    axis([0 360 0 90]); yticks(0:30:90); xticks(0:60:360);
elseif isnumeric(skyline)
    imagesc(flipud(1-skyline.'))
    colormap(gray);
    [n,m] = size(skyline);
    axis([0 n 0 m]);
    xticks(linspace(0,n,7)); xticklabels(linspace(0,360,7))
    yticks(linspace(0,m,4)); yticklabels(linspace(90,0,4))
end
font_size = 16;
xlabel(['Azimuth [' char(176) ']'], 'fontsize', font_size)
ylabel(['Altitude [' char(176) ']'], 'fontsize', font_size)
a = {'0 (N)';'60';'120';'180 (S)';'240';'300';'360 (N)'};
set(gca,'XTickLabel',a,'FontName','Calibri','fontsize',font_size)

end


function shading_factor = compute_shading_factor(skyline, sun_azimuth, ...
    sun_altitude)
%COMPUTE_SHADING_FACTOR Compute the shading factor depending on the sun
%position and the skyline of the PV system
%
% A linear sun trajectory is assumed for the calculations
% The skyline can be given as a cell or as a matrix

shading_factor = ones(length(sun_azimuth),1);
if iscell(skyline)   
    if all(size(skyline) == [1,2])
        x_azim = skyline{1,1}'; y_alti = skyline{1,2}';
    elseif all(size(skyline) == [1,1])
        skl = skyline{1,1};
        x_azim = skl(:,1); y_alti = skl(:,2);
    end
    [~,az_idx_prev] = min(abs(sun_azimuth(1) - x_azim));
    for i = 2:length(sun_azimuth)
        [~,az_idx] = min(abs(sun_azimuth(i) - x_azim));
        % Azimuth index will only be greater than its previous value when
        % crossing the north. In those cases, it's safe to assume that there
        % will be no direct irradiance, so this component is reduced
        if az_idx_prev<az_idx
            n = az_idx-az_idx_prev+1;
            range_sun_alt = linspace(sun_altitude(i-1),sun_altitude(i),n)';
            range_skl_alt = y_alti(az_idx_prev:az_idx);
            shading_factor(i) = 1-sum(range_sun_alt<range_skl_alt)/n;
        end
        az_idx_prev = az_idx;
    end
else
    [n,m] = size(skyline);
    az_idx_prev = round(1+sun_azimuth(1)/(360/(n-1)));
    for i = 2:length(sun_azimuth)
        az_idx = round(1+sun_azimuth(i)/(360/(n-1)));
        if az_idx_prev<az_idx
            az_idx_range = az_idx_prev:az_idx;
            len = length(az_idx_range);
            sun_alt_range = linspace(sun_altitude(i-1),sun_altitude(i),len);
            alt_idx_range = round(1+sun_alt_range/(90/(m-1)));
            SF = 0;
            for idx = [az_idx_range;alt_idx_range]
                SF = SF + skyline(idx(1),idx(2));
            end
            shading_factor(i) = 1 - SF/length(az_idx_range);
        end
        az_idx_prev = az_idx;
    end
end

end