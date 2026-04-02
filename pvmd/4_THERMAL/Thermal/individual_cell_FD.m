function cell_temp = individual_cell_FD(ambient_temp,wind_speed,cell_abs,...
    TOOLBOX_input,MODULE_output)
%INDIVIDUAL_CELL_FD Calculate cell temperature using fluid dynamic model
% 
% The cell temperature is calculated using the fluid-dynamic model as
% explained in the solar book Smets et al with improvements. The
% irradiance distribution of the cell is calculated by solving energy
% balance equations between three elements: front glass, solar cell, and
% rear glass.
%
% Parameters
% ----------
% ambient_temp : double
%   Ambient temperature at each time instant (celsius)
% wind_speed : double
%   Wind speed at each time instant (m/s)
% cell_abs : double
%   Absorbed irradiance by each cell at each time instant
% TOOLBOX_input : struct
%   Simulation parameters
% MODULE_output : struct
%   Results from the MODULE module
%
% Returns
% -------
% cell_temp : double
%   Cell temperature at each time instant (celsius)
%
% Performed by unknown (A. Jamodkar? E. Garcia?). Improved by A. Nour.
% Commented by A. Alcaniz

%% Load data

[pvmd_folder,~,~] = get_folder_structure;
file_path = fullfile(pvmd_folder,'4_THERMAL','data','thermal_params.mat');
load(file_path, 'anenometer_height','glass_conduct',...
    'glass_emissivity','boltzmann_ct')

%mounting height of the module convert [cm]->[m]
ModMountHeight=0.01*TOOLBOX_input.Scene.module_mounting.ModMountHeight;

% Efficiency of the PV module [-]
cell_eff = TOOLBOX_input.thermal.cell_eff;

% Temperature coefficient of the PV module [1/K]
temp_coeff=TOOLBOX_input.thermal.temp_coeff;

% Thickness of the glass [m]
glass_thickness=TOOLBOX_input.thermal.glass_thickness;

% number_cells      number of cells in a module
number_cells=MODULE_output.N;

% cell_area         area of the PV cell (cm2), convert from m2
cell_area=MODULE_output.A*10000;

% glass_area        area of the glass (cm2) per cell
glass_area=MODULE_output.Amod*10000/number_cells;

% tilt              tilt of the module (degrees)
tilt=MODULE_output.ModTilt;

%% Precalculate certain parameters to reduce calculations inside the loop
wind_speed_module = wind_speed*(ModMountHeight/anenometer_height)^0.2;
number_hours = size(ambient_temp,1);

% Absorbed irradiance for the front and rear glass
front_glass_abs = 0.014*cell_abs;
rear_glass_abs = 0.001*cell_abs;

% Thermal heat transfer
H_cond = (1/(1/(glass_conduct/glass_thickness)+1/(0.24/0.002))+(0.1578/0.002))*cell_area;

% Params for radiative heat transfer
F_sky_front = 0.5*(1 + cosd(tilt));
F_gro_front = 0.5*(1 - cosd(tilt));
F_sky_back = 0.5*(1 + cosd(180-tilt));
F_gro_back = 0.5*(1 - cosd(180-tilt));
radiative_ct = glass_emissivity*boltzmann_ct*glass_area;

%% Iterative loop for each hour for each cell in the module
% Initialize parameters to store the results
cell_temp = zeros(size(cell_abs));

for h = 1:number_hours
    ambient_temp_h = ambient_temp(h)+273.15;
    ground_temp_h = ambient_temp_h;
    sky_temp_h = 0.68*(0.0552*ambient_temp_h^1.5) + 0.32*ambient_temp_h;
    
    % Initial guess of temperatures of the cell and glasses
    cell_temp_h_guess = ambient_temp_h + 30;
    front_glass_temp = cell_temp_h_guess - 2;
    rear_glass_temp = cell_temp_h_guess - 2;

    for cell = 1:number_cells
        cell_temp_h = 0;          
        
        while abs(cell_temp_h_guess - cell_temp_h) >= 0.0001
            cell_temp_h = cell_temp_h_guess;
            
            % Heat Flux 
            Heat_flux = cell_abs(h,cell)*(1 - cell_eff*(1 + temp_coeff*(cell_temp_h - 300)));
            
            % Free convection heat transfer
            h_conv_free_front = 1.31*abs(front_glass_temp - ambient_temp_h)^(1/3);
            h_conv_free_back = 1.31*abs(rear_glass_temp - ambient_temp_h)^(1/3);
            h_conv_forced = 2.8 + 3*wind_speed_module(h);
            
            % Total convection
            H_conv_front = (h_conv_free_front^3 + h_conv_forced^3)^(1/3)*glass_area;
            H_conv_rear =  (h_conv_free_back^3 + h_conv_forced^3)^(1/3)*glass_area;
            
            % Radiative heat transfer
            h_rad_sky_F = radiative_ct*F_sky_front*(front_glass_temp + sky_temp_h)*(front_glass_temp^2 + sky_temp_h^2);
            h_rad_gro_F = radiative_ct*F_gro_front*(front_glass_temp + ground_temp_h)*(front_glass_temp^2 + ground_temp_h^2);
            h_rad_sky_B = radiative_ct*F_sky_back*(rear_glass_temp + sky_temp_h)*(rear_glass_temp^2 + sky_temp_h^2);
            h_rad_gro_B = radiative_ct*F_gro_back*(rear_glass_temp + ground_temp_h)*(rear_glass_temp^2 + ground_temp_h^2);

            % Front glass and rear glass temperature
            front_glass_temp = (ambient_temp_h*H_conv_front + sky_temp_h*h_rad_sky_F + ground_temp_h*h_rad_gro_F + cell_temp_h*H_cond + front_glass_abs(h,cell)*glass_area) / (H_conv_front + h_rad_sky_F + h_rad_gro_F + H_cond);
            rear_glass_temp = (ambient_temp_h*H_conv_rear + sky_temp_h*h_rad_sky_B + ground_temp_h*h_rad_gro_B + cell_temp_h*H_cond + rear_glass_abs(h,cell)*glass_area) / (H_conv_rear + h_rad_sky_B + h_rad_gro_B + H_cond);
            cell_temp_h_guess = (front_glass_temp*H_cond + rear_glass_temp*H_cond + cell_area*Heat_flux) / (2.0*H_cond);
        end
        cell_temp(h,cell) = cell_temp_h_guess - 273.15;
    end
    
end

end
