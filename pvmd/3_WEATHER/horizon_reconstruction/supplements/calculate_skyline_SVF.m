function [SVF, skyline] = calculate_skyline_SVF(h_iy,bbox_iy,pts,...
    azimuth,tilt,radius)
%CALCULATE_SKYLINE_SVF Compute the skyline and SVF
%
% Parameters
% ----------
% h_iy : double
%   Height map grid around the point of interest with z values in meters
% bbox_iy : double
%   2x2 matrix with the absolute borders of H in meters
% pts : double
%   Position [x,y,z] of the point of interest
% lat : double
%   Latitude of the point of interest
% lon : double
%   Longitude of the point of interest
% azimuth : double
%   Azimuth angle in degrees from north of the plane of array (POA)
% tilt : double
%   Tilt angle in degrees from flat of the plane of array (POA)
% radius : double
%   Perimeter around the point of interest to evaluate
%
% Returns
% -------
% SVF : double
%   Sky view factor of the point of interest
% skyline : cell
%   Skyline surrounding the point of interest

alt_steps = 90; % # of steps of sky dome in 90 degrees
azim_steps = 720; % # of steps of sky dome in 360 degrees

% Get the horizon elevation
skyline = horizon_scanner(h_iy,bbox_iy,pts,azimuth,tilt,radius,azim_steps);

% Generate sky sector view matrix
sky_vf = generate_sky_dome_matrix(alt_steps,azim_steps,azimuth,tilt);

% Get sky view factor
[SVF, skyline] = SVF_calculation(skyline, sky_vf);

end


function shapes = horizon_scanner(H,bbox,pts,poa_azim,poa_tilt,radius,...
    azim_steps)
%HORIZON_SCANNER scan the horizon to create the skyline
%
% The first version of Horizonscanner was written and used by Martijn
% Keijzer in his MSc Thesis about PV yield calculations in complex urban
% geometries. The same file was used by Josu Etxebarria in the rooftop PV
% potential tool. In this version the code has been rewritten by Maarten 
% Verkou to improve computation speed for the rooftop analysis tool.
%
% Parameters
% ----------
% H : double
%   Height map grid with z values in meters
% bbox : double
%   2x2 matrix with the absolute borders of H in meters
% pts : double
%   Position [x,y,z] of the point of interest
% poa_azim : double
%   Azimuth angle in degrees from north of the plane of array (POA)
% poa_tilt : double
%   Tilt angle in degrees from flat of the plane of array (POA)
% radius : double
%   Perimeter around pts to evaluate
% azimSteps : double
%   Number of azimuth steps for the 360 degrees view
%
% Returns
% -------
% skyline : cell
%   Skyline surrounding the point of interest

% Create a vector with the azimuth values
azimInt = 360/azim_steps;
azimVec = azimInt*((1:azim_steps)'-1.0);
shape_azim = azimVec;
max_r = 1/tand(3*360/azim_steps)*(2/1.75); % multislices

% Convert point matrix to single column vectors 
xlin = linspace(bbox(1,1),bbox(2,1),size(H,2)); % east->west (Low->High)
ylin = linspace(bbox(2,2),bbox(1,2),size(H,1)); % north->south (High->Low)
X = repmat(xlin,size(H,1),1);  x=X(:); % [x1 x1 x1 x2 x2 x2 ... xn xn xn]
y = repmat(ylin',size(H,2),1);  % [y1 y2 ym y1 y2 ym ... y1 y2 ym]
h = H(:); % [h11 h21 h31 hm1 h12 h22 h32 hm2 ... h1n h2n h3n hmn] 
xyz = [x,y,h]; 

% Calculate distances
dxyz = xyz - pts;
rho = sqrt((xyz(:,1)-pts(1)).^2+(xyz(:,2)-pts(2)).^2);
% Get points within the radius and higher than the module height
inrad = rho<radius & h>=pts(3);
xyz_inr = xyz(inrad,:); dxyz_inr = dxyz(inrad,:); rho_inr = rho(inrad);

% Remove points behind the plane of array, excluding nearby points
out = rho_inr>0.9*radius;
xyz_far = xyz(out,:);
A_far = mod(atan2d(dxyz_inr(out,1),dxyz_inr(out,2)),360);  
a_far = atand(dxyz_inr(out,3)./rho_inr(out));  
aoi = sind(poa_tilt).*sind(90-a_far).*cosd(poa_azim-A_far)+cosd(poa_tilt).*cosd(90-a_far)>=0;
xyz_poa = xyz_far(aoi,:); a_poa = a_far(aoi);
xyz_near = xyz_inr(~out,:);

% Remove far away points with altitude below 50 percentile 
a_perc = prctile(a_poa,50);
xyz_perc = xyz_poa(a_poa>a_perc,:);

xyz_tot = [xyz_perc;xyz_near];

shapes = cell(1);
horizon = zeros([azim_steps 1]);  % Create an emptyhorizon
% Get the relative cartesian coordinates with current point as origin
dxyh = xyz_tot(xyz_tot(:,3)>pts(3),:)-pts;

% Convert to spherical coordinates
rho_p = sqrt(dxyh(:,1).^2+dxyh(:,2).^2); 
if poa_tilt > 80
    idx = rho_p>=2.5;
    rho_p = rho_p(idx); dxyh = dxyh(idx,:);
end
if isempty(dxyh)
    dxyh = [0,10,0;10,0,0;0,-10,0;-10,0,0];
    rho_p = [10;10;10;10];
    disp('no points')
end
% Remove points very close to wall;
elev_p = atand(dxyh(:,3)./rho_p);  
azim_p = mod(atan2d(dxyh(:,1),dxyh(:,2)),360);  
% Round azimuth values to nearest points in azim vector
azim_idx = interp1(azimVec,find(azimVec==azimVec),azim_p,'nearest','extrap');
% Get a sorted matrix with corresponding elevation and distance
Tbl_sort = sortrows([azim_idx,elev_p,rho_p],1);
% Get the row indexes of certain azimuths
rowid = [0;find(diff(Tbl_sort(:,1)));size(Tbl_sort(:,1),1)];  

% Loop through the different azimuth angles
for k = 1:(numel(rowid)-1)
    f = (rowid(k)+1):(rowid(k+1)); % Index of current azimuth in azimVec
    [El,idx] = max(Tbl_sort(f,2));
    Rh = Tbl_sort(f(idx),3);
    Az = Tbl_sort(f(1),1); curEl = horizon(Az);
    horizon(Az) = max(El,curEl);
    if Rh < max_r
        % amount of slices to each side. 1.75 is middle of 2 & sqrt(2)
        angwidth = round((atand(1/Rh)/(360/azim_steps))/1.75);                                    
        % Make elevation the maximum within angwidth
        i1 = mod(Az-(angwidth+1),azim_steps)+1;
        i2 = mod(Az+(angwidth-1),azim_steps)+1;
        if i2<i1 % This happens when i1 < 0 (e.g. <360)
            horizon(i1:azim_steps) = max(El,horizon(i1:azim_steps));
            horizon(1:i2) = max(El,horizon(1:i2));
        else
            horizon(i1:i2) = max(El,horizon(i1:i2));
        end
    end
end

% Create the shape
shape_elev = horizon(:);
shapes{1} = [shape_azim shape_elev];

end


function sky_sector_vf = generate_sky_dome_matrix(rows,cols,poa_azim,poa_tilt)
%GENERATE_SKY_DOME_MATRIX Generate a matrix with the sky dome
%
% Parameters
% ----------
% rows : double
%   Number of altitude steps for the 360 degrees view
% cols : double
%   Number of azimuth steps for the 360 degrees view
% poa_azim : double
%   Azimuth angle in degrees from north of the plane of array (POA)
% poa_tilt : double
%   Tilt angle in degrees from flat of the plane of array (POA)
%
% Returns
% -------
% sky_sector_vf : double
%   Matrix with the view factor of the sky sector

% Create vectors of polar coordinates of the centre of dome elements
azimV = 360/cols*((1:cols)-0.5);
altV = 90/rows*((1:rows)'-0.5);
z = 90 - altV;

% Cosine of angle of incidence 
f_aoi = sind(poa_tilt)*sind(z)*cosd(poa_azim-azimV)+cosd(poa_tilt)*cosd(z);
f_aoi(f_aoi<0) = 0;

% Sine of zenith accounts for sky element size in dome
jacob = repmat(sind(z),1,cols);

% Sky sector view factor matrix following A.Calcabrini
cos_0 = repmat(cosd(z),1,cols);
sky_tot_max = sum(jacob.*cos_0,'all');
sky_sector_vf = (jacob.*f_aoi)/sky_tot_max;

end

function [SVF, skyline] = SVF_calculation(skyline_raw,sky_sector_vf)
% SVF_CALCULATION Calculate the sky view factor (SVF)
%
% Parameters
% ----------
% skyline_raw : cell
%   Skyline surrounding the point of interest
% sky_sector_vf : double
%   Matrix with the view factor of the sky sector
%
% Returns
% -------
% SVF : double
%   Sky view factor (SVF) of the point of interest
% skyline : cell
%   Skyline surrounding the point of interest

[ROWS,COLS] = size(sky_sector_vf);
azimV = 360/COLS*((1:COLS)-0.5);
altV = 90/ROWS*((1:ROWS)-0.5);

skyline_matrix = true(ROWS,COLS);
if iscell(skyline_raw)
    shapes_num = length(skyline_raw);
    if shapes_num>0
        for s=1:shapes_num
            skyline_matrix = skyline_matrix & ~inpolygon(...
                repmat(azimV,length(altV),1),...
                repmat(altV',1,length(azimV)),...
                [0;skyline_raw{s}(:,1);360],...
                [0;skyline_raw{s}(:,2);0]);
        end
    end
elseif size(skyline_raw,1) == COLS+1 && size(skyline_raw,2) == 2
    for iAz = 1:1:COLS
        idx1 = iAz;     % == find(skyline_raw(:,1)<azimV(iAz),1,'last'); 
        idx2 = iAz+1;   % == find(skyline_raw(:,1)>azimV(iAz),1,'first'); 
        skyline_matrix(:,iAz) = altV>((skyline_raw(idx1,2)+skyline_raw(idx2,2))/2);
    end
elseif size(skyline_raw,1) == ROWS && size(skyline_raw,2) == COLS
    skyline_matrix = skyline_matrix & skyline_raw;
else
    error('Invalid skyline!');
end

% Perform SVF calculations
SVF = sum(skyline_matrix.*sky_sector_vf,'all');

x_azim = linspace(0, 360, size(skyline_matrix,2));
y_alti = 90 - sum(skyline_matrix)/size(skyline_matrix,1)*90;
skyline = {[x_azim', y_alti']};
end
