function horizon_reconstruction(TOOLBOX_input)
%HORIZON_RECONSTRUCTION Reconstruct the horizon using LiDAR data
% 
% Given the position of a PV module, read LiDAR data around it and
% reconstruct the horizon in order to obtain the sky view factor (SVF) and
% the skyline
% 
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
%
% Based on RoofPotentialCalculation by M. Verkou. Adapted by A. Alcaniz for
% the toolbox

%---- Find the input dir with tiff files
[~,~,data_folder] = get_folder_structure;
input_dir = fullfile(data_folder,'Weather','Reconstructed Horizons',...
    'Height Files');

%---- Load characteristics of the point to be studied
latitude = TOOLBOX_input.irradiation.latitude;
longitude = TOOLBOX_input.irradiation.longitude;
radius = TOOLBOX_input.irradiation.radius;
if TOOLBOX_input.runPeriodic
    module_height = TOOLBOX_input.Scene.module_mounting.ModMountHeight/100; %from cm to m

    % The point of reference of the module azimuth is changed (from S=0, N=180)
    % to match the reference of the skyline (S=180, N=0)
    azimuth = TOOLBOX_input.Scene.module_mounting.ModAzimuth + 180;
    if azimuth >= 360; azimuth = azimuth - 360; end

else
    module_height = TOOLBOX_input.Scene.module_mounting(1).zCoordinate/100; %from cm to m

    % The point of reference of the module azimuth is changed (from S=0, N=180)
    % to match the reference of the skyline (S=180, N=0)
    azimuth = TOOLBOX_input.Scene.module_mounting(1).ModAzimuth + 180;
    if azimuth >= 360; azimuth = azimuth - 360; end
end

% The tilt is set to 0 because the effects of the tilt in the irradiance
% are already computed in the previous block
tilt = 0; 

% From WGS to Dutch National coords
[x, y] = wgs2rd(longitude, latitude);

%---- Import LiDAR 2D Raster
tiffs = dir([input_dir,filesep,'*.TIF']);
for i = 1:size(tiffs,1)
    [tiffs(i).H,~,~,tiffs(i).bbox] = geotiffread(fullfile(input_dir,tiffs(i).name));
end

%---- Merge/clip LiDAR data based on building contours
% Reduce points grid to points inside building bbox + radius + a bit extra
bbox_iy = [x,y;x,y] + 1.1*radius*[-1;1];
% Merge LiDAR data based on the box of interest
[h,bbox] = geotiffmerge(tiffs,bbox_iy);

%---- Obtain the height matrix of points surrounding the module
[h_iy,bbox_iy] = mapinpolygon(h,bbox,bbox_iy(:,1),bbox_iy(:,2),true);
h_iy(h_iy>10000) = NaN; % Remove noise (water reflect) points

%---- Obtain the height at the module location
bbox_mod = [x,y;x,y] + [-0.5;0.5];
[h_mod,~] = mapinpolygon(h,bbox,bbox_mod(:,1),bbox_mod(:,2),true);
% Height at the module location
h_location = h_mod(1,1);
% Considering the height input by the user relative to the ground
z = h_location + module_height;
pts = [x,y,z];

clear H bbox  % Clear the large files and bbox, to prevent further use

%---- SVF and skyline calculation

[SVF, skyline] = calculate_skyline_SVF(h_iy,bbox_iy,pts,...
    azimuth,tilt,radius);


%% Save the variables for later use
% The coordinates and city name are also saved for faster recognition of
% skyline
[~,~,data_folder] = get_folder_structure;
hor_rec_dir = fullfile(data_folder,'Weather','Reconstructed Horizons');
file_name = fullfile(hor_rec_dir,TOOLBOX_input.irradiation.skyline_file);
save(file_name,'skyline','SVF','latitude','longitude')

end
