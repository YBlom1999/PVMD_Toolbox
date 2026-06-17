function [pvt_collector_output] = ATES_storage_A(~,~, WEATHER_output)
%This function calculates the thermal & PV yield of water-based PVT collector integrated with ATES
% ATES
% Parameters
% ----------
% MODULE_output : struct
%   Output of this module block
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% pvt_collector_output : struct
%   Simulation results of the PVT module
% 
% Developed by A & ZUA.

I = mean(WEATHER_output.Irr,2)';
Tam = [WEATHER_output.ambient_temperature]';

idx = find(I~=0);                    % find the indices of non-zero elements in I
if length(idx) > length(I)           % if there are non-zero elements in I
    idx = idx(1:length(I));          % select the non-zero elements
end
I1 = I(idx);                         % select the corresponding elements from I

% Values for appartments
Appartments = 100;
Hospitals = 0;
Medium_office = 0;
Primary_school = 0;
Quickservice_restaurant = 0;
Restaurant = 0;
Secundairy_school = 0;
Small_hotel = 0;
Small_office = 0;
Standalone_retail = 0;
Warehouse = 0;
Large_hotel = 0;
Supermarket = 0;

% Values for PVT modules
PVT_Modules_per_appartment = 3;
PVT_Modules_per_Hospital = 50;
PVT_Modules_per_medium_office = 6;
PVT_Modules_per_primary_school = 8;
PVT_Modules_per_Quickservice_restaurant = 7;
PVT_Modules_per_Restaurant = 15;
PVT_Modules_per_Secundairy_school = 25;
PVT_Modules_per_Small_hotel = 12;
PVT_Modules_per_Small_office = 5;
PVT_Modules_per_Standalone_retail = 6;
PVT_Modules_per_Warehouse = 6;
PVT_Modules_per_Large_hotel = 40;
PVT_Modules_per_Supermarket = 25;

% Scaling to heat demand
heating_demand_per_appartment = 6000;  % kWh/year per house, around 6000 for average house
cooling_demand_per_appartment = 2100;  % kWh/year per house, around 2100 for average house
heating_demand_per_Hospital = 80000;
cooling_demand_per_Hospital = 20000;
heating_demand_per_medium_office = 20000;
cooling_demand_per_medium_office = 5000;
heating_demand_per_primary_school = 15000;
cooling_demand_per_primary_school = 6000;
heating_demand_per_Quickservice_restaurant = 3700;
cooling_demand_per_Quickservice_restaurant = 2000;
heating_demand_per_Restaurant = 7400;
cooling_demand_per_Restaurant = 4000;
heating_demand_per_Secundairy_school = 18200;
cooling_demand_per_Secundairy_school = 8000;
heating_demand_per_Small_hotel = 9000;
cooling_demand_per_Small_hotel = 10000;
heating_demand_per_Small_office = 11500;
cooling_demand_per_Small_office = 3000;
heating_demand_per_Standalone_retail = 3600;
cooling_demand_per_Standalone_retail = 2000;
heating_demand_per_Warehouse = 8500;
cooling_demand_per_Warehouse = 4000;
heating_demand_per_Large_hotel = 40000;
cooling_demand_per_Large_hotel = 10000;
heating_demand_per_Supermarket = 20000;
cooling_demand_per_Supermarket = 7500;

% Values for ATES system
Max_temp_diff = 14;          % Maximum fluctuation in temperature of aquifer
Min_temp_diff2 = 7;          % Minimum fluctuation in temperature of aquifer
% Min_temp_diff = 6;         % Minimum temperature difference between aquifers
Average_temp_diff = 1;       % Difference average temperature beginning and end year
Ending_temp_diff_1 = 3;      % Difference temperature aquifer 1 beginneng and end of year
Ending_temp_diff_2 = 3;      % Difference temperature aquifer 2 beginning and end of year
Max_temp = 45;               % Maximum temperature of the hottest aquifer
Min_temp = 5;                % Minimum temperature of the coldest aquifer

% Parameters
rho_w = 1000;                % Density (kg/m^3)
C_w = 4186;                  % Specific heat capacity (J/kg/K)
eta_hex = 0.4;

% Calculate the total heat demand
total_heat_demand = heating_demand_per_appartment * Appartments + ...
                    heating_demand_per_Hospital * Hospitals + ...
                    heating_demand_per_medium_office * Medium_office + ...
                    heating_demand_per_primary_school * Primary_school + ...
                    heating_demand_per_Quickservice_restaurant * Quickservice_restaurant + ...
                    heating_demand_per_Restaurant * Restaurant + ...
                    heating_demand_per_Secundairy_school * Secundairy_school + ...
                    heating_demand_per_Small_hotel * Small_hotel + ...
                    heating_demand_per_Small_office * Small_office + ...
                    heating_demand_per_Standalone_retail * Standalone_retail + ...
                    heating_demand_per_Warehouse * Warehouse + ...
                    heating_demand_per_Supermarket * Supermarket + ...
                    heating_demand_per_Large_hotel * Large_hotel; % kWh/year in total

% Calculate the total cooling demand
total_cooling_demand = cooling_demand_per_appartment * Appartments + ...
                       cooling_demand_per_Hospital * Hospitals + ...
                       cooling_demand_per_medium_office * Medium_office + ...
                       cooling_demand_per_primary_school * Primary_school + ...
                       cooling_demand_per_Quickservice_restaurant * Quickservice_restaurant + ...
                       cooling_demand_per_Restaurant * Restaurant + ...
                       cooling_demand_per_Secundairy_school * Secundairy_school + ...
                       cooling_demand_per_Small_hotel * Small_hotel + ...
                       cooling_demand_per_Small_office * Small_office + ...
                       cooling_demand_per_Standalone_retail * Standalone_retail + ...
                       cooling_demand_per_Warehouse * Warehouse + ...
                       cooling_demand_per_Supermarket * Supermarket + ...
                       cooling_demand_per_Large_hotel * Large_hotel; % kWh/year in total

% Parameters Aquifer in general 
% k_soil = 1.5;               % Thermal Conductivity of Soil (W/(m*K))
% T_soil= 12 + 273.15;        % °K
n_porosity = 0.15;            % Porosity aquifer 

% Parameters Aquifer 1 
T_initial_1 = 25+273.15;      % °K
V_1 = 1.5e4;
m_1 = rho_w * V_1;            % Mass of First Aquifer (kg)

% Parameters Aquifer 2 
T_initial_2 = 15+273.15;      % °K
V_2 = 1.25e4;
m_2 = rho_w * V_2;            % Mass of second Aquifer (kg)

% Timesteps
hours_year = linspace(0, length(Tam), length(Tam));

%% Theoretical functions for cool and heat demand

% Assuming 'hours_year' is defined as an array representing each hour of the year
cooling_demand_yearly = (total_cooling_demand/24) * sin(pi * (hours_year/24 - 90) / 180);
cooling_demand_yearly(cooling_demand_yearly < 0) = 0;
heating_demand_yearly = (total_heat_demand/24) * sin(pi * (hours_year/24 - 270) / 180);
heating_demand_yearly(heating_demand_yearly < 0) = 0;

% Normalize the heating demand
total_current_heating_demand = sum(heating_demand_yearly);
normalized_heating_demand_yearly = (total_heat_demand) / total_current_heating_demand * heating_demand_yearly;

% Normalize the cooling demand
total_current_cooling_demand = sum(cooling_demand_yearly);
normalized_cooling_demand_yearly = (total_cooling_demand) / total_current_cooling_demand * cooling_demand_yearly;

Total_yearly_heat_demand = sum(normalized_heating_demand_yearly);
Total_yearly_cold_demand = sum(normalized_cooling_demand_yearly);

% Create a sinusoidal function for the entire year, repeating every 24 hours
omega = 2 * pi / 24;        % daily frequency
amplitude = 0.2;            % Modify as needed for variability
shift = 1;                  % To ensure the average remains the same
daily_variation_yearly = amplitude * sin(omega * mod(hours_year, 24)) + shift;

% Apply the daily variation to the entire year's heating and cooling demand
normalized_heating_demand_yearly = normalized_heating_demand_yearly .* daily_variation_yearly;
normalized_cooling_demand_yearly = normalized_cooling_demand_yearly .* daily_variation_yearly;

%% Everything together
% Inputs
% Cooling season
Injected_heat_cooling_season = normalized_cooling_demand_yearly * (eta_hex);      
Extracted_cold_cooling_season = normalized_cooling_demand_yearly / eta_hex;               
% Heating season
Injected_cold_heating_season = normalized_heating_demand_yearly * (eta_hex);
Extracted_heat_heating_season = normalized_heating_demand_yearly / eta_hex;

% Not below zero
Injected_heat_cooling_season(Injected_heat_cooling_season < 0) = 0;
Extracted_cold_cooling_season(Extracted_cold_cooling_season < 0) = 0;
Injected_cold_heating_season(Injected_cold_heating_season < 0) = 0;
Extracted_heat_heating_season(Extracted_heat_heating_season < 0) = 0;

% Sum inputs
Total_injected_heat = sum(Injected_heat_cooling_season);
Total_extracted_cold = sum(Extracted_cold_cooling_season);
Total_injected_cold = sum(Injected_cold_heating_season);
Total_extracted_heat = sum(Extracted_heat_heating_season);

%% Actual heat demand
Heat_demand_midrise_appartment = readmatrix('Heat_demand_midrise_appartment.xlsx', 'Range', 'G2:G8761');     % in J')
Heat_demand_midrise_appartment_kwh = Heat_demand_midrise_appartment / 3.6e6; % in Kwh
% Sum_midrise_appartment = sum(Heat_demand_midrise_appartment_kwh);
total_current_midrise_appartment_demand = sum(Heat_demand_midrise_appartment_kwh);
normalized_Heat_demand_midrise_appartment = ((heating_demand_per_appartment / total_current_midrise_appartment_demand) * Heat_demand_midrise_appartment_kwh);
normalized_Heat_demand_midrise_appartment(isnan(normalized_Heat_demand_midrise_appartment)) = 0;
% Check_midrise_appartment = sum(normalized_Heat_demand_midrise_appartment);

Heat_demand_Hospital = readmatrix('Heat_demand_hospital.xlsx', 'Range', 'G2:G8761');
Heat_demand_Hospital_kwh = Heat_demand_Hospital / 3.6e6 ;
total_current_hospital_demand = sum(Heat_demand_Hospital_kwh);
normalized_Heat_demand_Hospital = (heating_demand_per_Hospital / total_current_hospital_demand) * Heat_demand_Hospital_kwh;
normalized_Heat_demand_Hospital(isnan(normalized_Heat_demand_Hospital)) = 0;
% Check_hospital = sum(normalized_Heat_demand_Hospital);

Heat_demand_Large_hotel = readmatrix('Heat_demand_Large_hotel.xlsx', 'Range', 'G2:G8761');
Heat_demand_Large_hotel_kwh = Heat_demand_Large_hotel / 3.6e6 ;
total_current_large_hotel_demand = sum(Heat_demand_Large_hotel_kwh);
normalized_Heat_demand_Large_hotel = (heating_demand_per_Large_hotel / total_current_large_hotel_demand) * Heat_demand_Large_hotel_kwh;
normalized_Heat_demand_Large_hotel(isnan(normalized_Heat_demand_Large_hotel)) = 0;
% Check_large_hotel = sum(normalized_Heat_demand_Large_hotel);

Heat_demand_Medium_office = readmatrix('Heat_demand_Medium_office.xlsx', 'Range', 'G2:G8761');
Heat_demand_Medium_office_kwh = Heat_demand_Medium_office / 3.6e6;
total_current_medium_office_demand = sum(Heat_demand_Medium_office_kwh);
normalized_Heat_demand_Medium_office = (heating_demand_per_medium_office / total_current_medium_office_demand) * Heat_demand_Medium_office_kwh;
normalized_Heat_demand_Medium_office(isnan(normalized_Heat_demand_Medium_office)) = 0;
% Check_medium_office = sum(normalized_Heat_demand_Medium_office);

Heat_demand_Primary_school = readmatrix('Heat_demand_Primary_school.xlsx', 'Range', 'G2:G8761');
Heat_demand_Primary_school_kwh = Heat_demand_Primary_school / 3.6e6 ;
total_current_primary_school_demand = sum(Heat_demand_Primary_school_kwh);
normalized_Heat_demand_Primary_school = (heating_demand_per_primary_school / total_current_primary_school_demand) * Heat_demand_Primary_school_kwh;
normalized_Heat_demand_Primary_school(isnan(normalized_Heat_demand_Primary_school)) = 0;
% Check_primary_school = sum(normalized_Heat_demand_Primary_school);

Heat_demand_Quickservicerestaurant = readmatrix('Heat_demand_QuickServiceRestaurant.xlsx', 'Range', 'G2:G8761');
Heat_demand_Quickservicerestaurant_kwh = Heat_demand_Quickservicerestaurant / 3.6e6 ;
total_current_Quickservicerestaurant_demand = sum(Heat_demand_Quickservicerestaurant_kwh);
normalized_Heat_demand_Quickservicerestaurant = (heating_demand_per_Quickservice_restaurant / total_current_Quickservicerestaurant_demand) * Heat_demand_Quickservicerestaurant_kwh;
normalized_Heat_demand_Quickservicerestaurant(isnan(normalized_Heat_demand_Quickservicerestaurant)) = 0;
% Check_Quickservicerestaurant = sum(normalized_Heat_demand_Quickservicerestaurant);

Heat_demand_Restaurant = readmatrix('Heat_demand_restaurant.xlsx', 'Range', 'G2:G8761');
Heat_demand_Restaurant_kwh = Heat_demand_Restaurant / 3.6e6 ;
total_current_restaurant_demand = sum(Heat_demand_Restaurant_kwh);
normalized_Heat_demand_Restaurant = (heating_demand_per_Restaurant / total_current_restaurant_demand) * Heat_demand_Restaurant_kwh;
normalized_Heat_demand_Restaurant(isnan(normalized_Heat_demand_Restaurant)) = 0;
% Check_restaurant = sum(normalized_Heat_demand_Restaurant);

Heat_demand_Secundairy_school = readmatrix('Heat_demand_secundairy_school.xlsx', 'Range', 'G2:G8761');
Heat_demand_Secundairy_school_kwh = Heat_demand_Secundairy_school / 3.6e6 ;
total_current_secundairy_school_demand = sum(Heat_demand_Secundairy_school_kwh);
normalized_Heat_demand_Secundairy_school = (heating_demand_per_Secundairy_school / total_current_secundairy_school_demand) * Heat_demand_Secundairy_school_kwh;
normalized_Heat_demand_Secundairy_school(isnan(normalized_Heat_demand_Secundairy_school)) = 0;
% Check_secundairy_school = sum(normalized_Heat_demand_Secundairy_school);

Heat_demand_Small_hotel = readmatrix('Heat_demand_small_hotel.xlsx', 'Range', 'G2:G8761');
Heat_demand_Small_hotel_kwh = Heat_demand_Small_hotel / 3.6e6 ;
total_current_small_hotel_demand = sum(Heat_demand_Small_hotel_kwh);
normalized_Heat_demand_Small_hotel = (heating_demand_per_Small_hotel / total_current_small_hotel_demand) * Heat_demand_Small_hotel_kwh;
normalized_Heat_demand_Small_hotel(isnan(normalized_Heat_demand_Small_hotel)) = 0;
% Check_small_hotel = sum(normalized_Heat_demand_Small_hotel);

Heat_demand_Small_office = readmatrix('Heat_demand_small_office.xlsx', 'Range', 'G2:G8761');
Heat_demand_Small_office_kwh = Heat_demand_Small_office / 3.6e6 ;
total_current_small_office_demand = sum(Heat_demand_Small_office_kwh);
normalized_Heat_demand_Small_office = (heating_demand_per_Small_office / total_current_small_office_demand) * Heat_demand_Small_office_kwh;
normalized_Heat_demand_Small_office(isnan(normalized_Heat_demand_Small_office)) = 0;
% Check_small_office = sum(normalized_Heat_demand_Small_office);

Heat_demand_Standalone_retail = readmatrix('Heat_demand_standalone_retail.xlsx', 'Range', 'G2:G8761');
Heat_demand_Standalone_retail_kwh = Heat_demand_Standalone_retail / 3.6e6 ;
total_current_standalone_retail_demand = sum(Heat_demand_Standalone_retail_kwh);
normalized_Heat_demand_Standalone_retail = (heating_demand_per_Standalone_retail / total_current_standalone_retail_demand) * Heat_demand_Standalone_retail_kwh;
normalized_Heat_demand_Standalone_retail(isnan(normalized_Heat_demand_Standalone_retail)) = 0;
% Check_standalone_retail = sum(normalized_Heat_demand_Standalone_retail);

Heat_demand_Warehouse = readmatrix('Heat_demand_warehouse.xlsx', 'Range', 'G2:G8761');
Heat_demand_Warehouse_kwh = Heat_demand_Warehouse / 3.6e6 ;
total_current_warehouse_demand = sum(Heat_demand_Warehouse_kwh);
normalized_Heat_demand_Warehouse = (heating_demand_per_Warehouse / total_current_warehouse_demand) * Heat_demand_Warehouse_kwh;
normalized_Heat_demand_Warehouse(isnan(normalized_Heat_demand_Warehouse)) = 0;
% Check_warehouse = sum(normalized_Heat_demand_Warehouse);

Heat_demand_Supermarket = readmatrix('Heat_demand_supermarket1.xlsx', 'Range', 'G2:G8761');     % in J')
Heat_demand_Supermarket_kwh = Heat_demand_Supermarket / 3.6e6; % in Kwh
% Sum_supermarket = sum(Heat_demand_Supermarket_kwh);
total_current_Supermarket_demand = sum(Heat_demand_Supermarket_kwh);
normalized_Heat_demand_Supermarket = ((heating_demand_per_Supermarket / total_current_Supermarket_demand) * Heat_demand_Supermarket_kwh);
normalized_Heat_demand_Supermarket(isnan(normalized_Heat_demand_Supermarket)) = 0;
% Check_Supermarket = sum(normalized_Heat_demand_Supermarket);

% Summing the normalized heat demands of all building types
Total_heat_demand_all_buildings = normalized_Heat_demand_midrise_appartment * Appartments + ...
                                  normalized_Heat_demand_Hospital * Hospitals+ ...
                                  normalized_Heat_demand_Large_hotel * Large_hotel + ...
                                  normalized_Heat_demand_Medium_office * Medium_office+ ...
                                  normalized_Heat_demand_Primary_school * Primary_school+ ...
                                  normalized_Heat_demand_Quickservicerestaurant * Quickservice_restaurant+ ...
                                  normalized_Heat_demand_Restaurant * Restaurant+ ...
                                  normalized_Heat_demand_Secundairy_school * Secundairy_school+ ...
                                  normalized_Heat_demand_Small_hotel * Small_hotel+ ...
                                  normalized_Heat_demand_Small_office * Small_office+ ...
                                  normalized_Heat_demand_Standalone_retail * Standalone_retail+ ...
                                  normalized_Heat_demand_Supermarket * Supermarket+ ...
                                  normalized_Heat_demand_Warehouse * Warehouse;
Total_heat_demand_all_buildings = Total_heat_demand_all_buildings';
Total_annual_heat_demand = sum(Total_heat_demand_all_buildings);

%% Actual cold demand
Cold_demand_midrise_appartment = readmatrix('Heat_demand_midrise_appartment.xlsx', 'Range', 'H2:H8761');     % in J')
Cold_demand_midrise_appartment_kwh = Cold_demand_midrise_appartment / 3.6e6; % in Kwh
% Sum_midrise_appartment_Cold = sum(Cold_demand_midrise_appartment_kwh);
total_current_midrise_appartment_demand = sum(Cold_demand_midrise_appartment_kwh);
normalized_Cold_demand_midrise_appartment = ((cooling_demand_per_appartment / total_current_midrise_appartment_demand) * Cold_demand_midrise_appartment_kwh);
normalized_Cold_demand_midrise_appartment(isnan(normalized_Cold_demand_midrise_appartment)) = 0;
% Check_midrise_appartment_Cold = sum(normalized_Cold_demand_midrise_appartment);

Cold_demand_Hospital = readmatrix('Heat_demand_hospital.xlsx', 'Range', 'H2:H8761');
Cold_demand_Hospital_kwh = Cold_demand_Hospital / 3.6e6 ;
total_current_hospital_demand_Cold = sum(Cold_demand_Hospital_kwh);
normalized_Cold_demand_Hospital = (cooling_demand_per_Hospital / total_current_hospital_demand_Cold) * Cold_demand_Hospital_kwh;
normalized_Cold_demand_Hospital(isnan(normalized_Cold_demand_Hospital)) = 0;
% Check_hospital_cold = sum(normalized_Cold_demand_Hospital);

Cold_demand_Large_hotel = readmatrix('Heat_demand_Large_hotel.xlsx', 'Range', 'H2:H8761');
Cold_demand_Large_hotel_kwh = Cold_demand_Large_hotel / 3.6e6 ;
total_current_large_hotel_demand = sum(Cold_demand_Large_hotel_kwh);
normalized_Cold_demand_Large_hotel = (cooling_demand_per_Large_hotel / total_current_large_hotel_demand) * Cold_demand_Large_hotel_kwh;
normalized_Cold_demand_Large_hotel(isnan(normalized_Cold_demand_Large_hotel)) = 0;
% Check_large_hotel_Cold = sum(normalized_Cold_demand_Large_hotel);

Cold_demand_Medium_office = readmatrix('Heat_demand_Medium_office.xlsx', 'Range', 'H2:H8761');
Cold_demand_Medium_office_kwh = Cold_demand_Medium_office / 3.6e6;
total_current_medium_office_demand = sum(Cold_demand_Medium_office_kwh);
normalized_Cold_demand_Medium_office = (cooling_demand_per_medium_office / total_current_medium_office_demand) * Cold_demand_Medium_office_kwh;
normalized_Cold_demand_Medium_office(isnan(normalized_Cold_demand_Medium_office)) = 0;
% Check_medium_office_Cold = sum(normalized_Cold_demand_Medium_office);

Cold_demand_Primary_school = readmatrix('Heat_demand_Primary_school.xlsx', 'Range', 'H2:H8761');
Cold_demand_Primary_school_kwh = Cold_demand_Primary_school / 3.6e6 ;
total_current_primary_school_demand = sum(Cold_demand_Primary_school_kwh);
normalized_Cold_demand_Primary_school = (cooling_demand_per_primary_school / total_current_primary_school_demand) * Cold_demand_Primary_school_kwh;
normalized_Cold_demand_Primary_school(isnan(normalized_Cold_demand_Primary_school)) = 0;
% Check_primary_school_Cold = sum(normalized_Cold_demand_Primary_school);

Cold_demand_Quickservicerestaurant = readmatrix('Heat_demand_QuickServiceRestaurant.xlsx', 'Range', 'H2:H8761');
Cold_demand_Quickservicerestaurant_kwh = Cold_demand_Quickservicerestaurant / 3.6e6 ;
total_current_Quickservicerestaurant_demand = sum(Cold_demand_Quickservicerestaurant_kwh);
normalized_Cold_demand_Quickservicerestaurant = (cooling_demand_per_Quickservice_restaurant / total_current_Quickservicerestaurant_demand) * Cold_demand_Quickservicerestaurant_kwh;
normalized_Cold_demand_Quickservicerestaurant(isnan(normalized_Cold_demand_Quickservicerestaurant)) = 0;
% Check_Quickservicerestaurant_Cold = sum(normalized_Cold_demand_Quickservicerestaurant);

Cold_demand_Restaurant = readmatrix('Heat_demand_restaurant.xlsx', 'Range', 'H2:H8761');
Cold_demand_Restaurant_kwh = Cold_demand_Restaurant / 3.6e6 ;
total_current_restaurant_demand = sum(Cold_demand_Restaurant_kwh);
normalized_Cold_demand_Restaurant = (cooling_demand_per_Restaurant / total_current_restaurant_demand) * Cold_demand_Restaurant_kwh;
normalized_Cold_demand_Restaurant(isnan(normalized_Cold_demand_Restaurant)) = 0;
% Check_restaurant_Cold = sum(normalized_Cold_demand_Restaurant);

Cold_demand_Secundairy_school = readmatrix('Heat_demand_secundairy_school.xlsx', 'Range', 'H2:H8761');
Cold_demand_Secundairy_school_kwh = Cold_demand_Secundairy_school / 3.6e6 ;
total_current_secundairy_school_demand = sum(Cold_demand_Secundairy_school_kwh);
normalized_Cold_demand_Secundairy_school = (cooling_demand_per_Secundairy_school / total_current_secundairy_school_demand) * Cold_demand_Secundairy_school_kwh;
normalized_Cold_demand_Secundairy_school(isnan(normalized_Cold_demand_Secundairy_school)) = 0;
% Check_secundairy_school_Cold = sum(normalized_Cold_demand_Secundairy_school);

Cold_demand_Small_hotel = readmatrix('Heat_demand_small_hotel.xlsx', 'Range', 'H2:H8761');
Cold_demand_Small_hotel_kwh = Cold_demand_Small_hotel / 3.6e6 ;
total_current_small_hotel_demand = sum(Cold_demand_Small_hotel_kwh);
normalized_Cold_demand_Small_hotel = (cooling_demand_per_Small_hotel / total_current_small_hotel_demand) * Cold_demand_Small_hotel_kwh;
normalized_Cold_demand_Small_hotel(isnan(normalized_Cold_demand_Small_hotel)) = 0;
% Check_small_hotel_Cold = sum(normalized_Cold_demand_Small_hotel);

Cold_demand_Small_office = readmatrix('Heat_demand_small_office.xlsx', 'Range', 'H2:H8761');
Cold_demand_Small_office_kwh = Cold_demand_Small_office / 3.6e6 ;
total_current_small_office_demand = sum(Cold_demand_Small_office_kwh);
normalized_Cold_demand_Small_office = (cooling_demand_per_Small_office / total_current_small_office_demand) * Cold_demand_Small_office_kwh;
normalized_Cold_demand_Small_office(isnan(normalized_Cold_demand_Small_office)) = 0;
% Check_small_office_Cold = sum(normalized_Cold_demand_Small_office);

Cold_demand_Standalone_retail = readmatrix('Heat_demand_standalone_retail.xlsx', 'Range', 'H2:H8761');
Cold_demand_Standalone_retail_kwh = Cold_demand_Standalone_retail / 3.6e6 ;
total_current_standalone_retail_demand = sum(Cold_demand_Standalone_retail_kwh);
normalized_Cold_demand_Standalone_retail = (cooling_demand_per_Standalone_retail / total_current_standalone_retail_demand) * Cold_demand_Standalone_retail_kwh;
normalized_Cold_demand_Standalone_retail(isnan(normalized_Cold_demand_Standalone_retail)) = 0;
% Check_standalone_retail_Cold = sum(normalized_Cold_demand_Standalone_retail);

Cold_demand_Warehouse = readmatrix('Heat_demand_warehouse.xlsx', 'Range', 'H2:H8761');
Cold_demand_Warehouse_kwh = Cold_demand_Warehouse / 3.6e6 ;
total_current_warehouse_demand = sum(Cold_demand_Warehouse_kwh);
normalized_Cold_demand_Warehouse = (cooling_demand_per_Warehouse / total_current_warehouse_demand) * Cold_demand_Warehouse_kwh;
normalized_Cold_demand_Warehouse(isnan(normalized_Cold_demand_Warehouse)) = 0;
% Check_warehouse_Cold = sum(normalized_Cold_demand_Warehouse);

Cold_demand_Supermarket = readmatrix('Heat_demand_supermarket1.xlsx', 'Range', 'H2:H8761');     % in J')
Cold_demand_Supermarket_kwh = Cold_demand_Supermarket / 3.6e6; % in Kwh
% Sum_supermarket = sum(Cold_demand_Supermarket_kwh);
total_current_Supermarket_demand = sum(Cold_demand_Supermarket_kwh);
normalized_Cold_demand_Supermarket = ((cooling_demand_per_Supermarket / total_current_Supermarket_demand) * Cold_demand_Supermarket_kwh);
normalized_Cold_demand_Supermarket(isnan(normalized_Cold_demand_Supermarket)) = 0;
% Check_Supermarket_Cold = sum(normalized_Cold_demand_Supermarket);

% Summing the normalized cold demands of all building types
Total_cold_demand_all_buildings = normalized_Cold_demand_midrise_appartment * Appartments+ ...
                                  normalized_Cold_demand_Hospital * Hospitals+ ...
                                  normalized_Cold_demand_Large_hotel * Large_hotel+ ...
                                  normalized_Cold_demand_Medium_office * Medium_office+ ...
                                  normalized_Cold_demand_Primary_school * Primary_school+ ...
                                  normalized_Cold_demand_Quickservicerestaurant * Quickservice_restaurant + ...
                                  normalized_Cold_demand_Restaurant * Restaurant+ ...
                                  normalized_Cold_demand_Secundairy_school * Secundairy_school+ ...
                                  normalized_Cold_demand_Small_hotel * Small_hotel+ ...
                                  normalized_Cold_demand_Small_office * Small_office+ ...
                                  normalized_Cold_demand_Standalone_retail * Standalone_retail+ ...
                                  normalized_Cold_demand_Supermarket * Supermarket+ ...
                                  normalized_Cold_demand_Warehouse * Warehouse;
Total_cold_demand_all_buildings = Total_cold_demand_all_buildings';

%% Actual generation
A_g = 5.1;                      % Gross Collector Area (m2)
dm_f = 0.027;                   % Fluid Mass Flowrate (kg/s)

% Second Order Thermal Efficiency Equation unique coefficients
a1a = 11.294;                   % unique coefficient W/(m2K)                                          
a2a = 44.25;                    % unique coefficient W/(m2K2) 
etao_tha = 0.6853;              % unique coefficient
% Electrical Efficiency Equation unique coefficients
e1a = 0.3269;                   % unique coefficient W/(m2K) 
etao_ele = 0.1286;              % unique coefficient

% Calculate total PVT modules
total_PVT_modules = Appartments * PVT_Modules_per_appartment + ...
                    Hospitals * PVT_Modules_per_Hospital + ...
                    Medium_office * PVT_Modules_per_medium_office + ...
                    Primary_school * PVT_Modules_per_primary_school + ...
                    Quickservice_restaurant * PVT_Modules_per_Quickservice_restaurant + ...
                    Restaurant * PVT_Modules_per_Restaurant + ...
                    Secundairy_school * PVT_Modules_per_Secundairy_school + ...
                    Small_hotel * PVT_Modules_per_Small_hotel + ...
                    Small_office * PVT_Modules_per_Small_office + ...
                    Standalone_retail * PVT_Modules_per_Standalone_retail + ...
                    Warehouse * PVT_Modules_per_Warehouse + ...
                    Supermarket * PVT_Modules_per_Supermarket + ...
                    Large_hotel * PVT_Modules_per_Large_hotel;

% Output total PVT modules
fprintf('Total PVT Modules: %d\n', total_PVT_modules);

%% Actual Reservoir temperature Evolution
% Preallocate temperature arrays
T_aquifer_1_actual = zeros(1, length(hours_year));
T_aquifer_2_actual = zeros(1, length(hours_year));

% Set initial temperatures
T_aquifer_1_actual(1) = T_initial_1;
T_aquifer_2_actual(1) = T_initial_2;

% Preallocate arrays for T_inflow and T_outflow
T_inflow_array = zeros(1, length(hours_year));
T_outflow_array = zeros(1, length(hours_year));
T_aquifer_1_celsius_actual = zeros(1, length(hours_year));
T_aquifer_2_celsius_actual = zeros(1, length(hours_year));
Qa_array = zeros(1, length(hours_year));
Qe_array = zeros(1, length(hours_year));
eta_E_array = zeros(1, length(hours_year));
eta_tha = zeros(1, length(hours_year));
tmp_array = zeros(1, length(hours_year));
% eta_E_new_array = zeros(1, length(hours_year));

% Set initial temperatures for the arrays
T_aquifer_1_celsius_actual(1) = T_initial_1 - 273.15; % Assuming T_initial_1 is in Kelvin
T_aquifer_2_celsius_actual(1) = T_initial_2 - 273.15; % Assuming T_initial_2 is in Kelvin

% Initial value for T_outflow
T_outflow = T_initial_1 - 273.15; 

% Calculate efficiencys aquifers
L_screen_1 = 1.02 * V_1^(1/3);
Rh1 = sqrt(V_1/(n_porosity*pi*L_screen_1));
Rth_1 = 0.66*Rh1;
A_1 = V_1 * (2/L_screen_1 + 2/Rth_1);
r_1 = (2 * V_1) / A_1;
h_1 = A_1 / (2 * pi * r_1);
x1 = A_1/V_1;
thermal_efficiency1 = 1-(-1.50*x1+0.6);

L_screen_2 = 1.02 * V_2^(1/3);
Rh2 = sqrt(V_2/(n_porosity*pi*L_screen_2));
Rth_2 = 0.66*Rh2;
A_2 = V_2 * (2/L_screen_2 + 2/Rth_2);
r_2 = (2 * V_2) / A_2;
h_2 = A_1 / (2 * pi * r_2);
x2 = A_2/V_2;
thermal_efficiency2 = 1-(-1.94*x2+0.6);

% Display the results for optimal radius and 
fprintf('Aquifer 1: Optimal Radius: %.2f meters, Optimal Height: %.2f meters\n', r_1, h_1);
fprintf('Aquifer 2: Optimal Radius: %.2f meters, Optimal Height: %.2f meters\n', r_2, h_2);

% Time loop
for i = 1:(length(hours_year))
    % Updat-e T_i for t8he current iteration
    if i == 1
        T_inflow = T_initial_2 - 273.15;  
        T_i = T_inflow;
        tmp = 0;
    else
        T_inflow = T_outflow - (T_outflow - T_aquifer_2_celsius_actual(i-1)) * (1-eta_hex);
        tmp = (T_i - T_aquifer_2_celsius_actual(i-1)) ./ I(i);
    end

    % Second Order Thermal Efficiency Equation
    eta_T = max(min(etao_tha - a1a * tmp - a2a * tmp.^2, 1), 0);
    % Electrical Efficiency Equation
    eta_E = max(min(etao_ele - e1a * tmp, 0.1286), 0);
    % Instantaneous thermal power at normal incidence
    Qa = eta_T .* A_g .* I(i);
    % Instantaneous electrical power at normal incidence
    Qe = eta_E .* A_g .* I(i);
    PVT_output_actual = size(I);
    % Update T_outflow based on T_inflow and Qa
    T_outflow = T_inflow + Qa / (C_w * dm_f);
    PVT_output_actual(i) = Qa * total_PVT_modules / 1000; % In Kw per hour for the entire year
    % Store T_inflow and T_outflow values in the arrays
    T_inflow_array(i) = T_inflow;
    T_outflow_array(i) = T_outflow;
    Qa_array(i) = Qa;
    Qe_array(i) = Qe;
    eta_E_array(i) = eta_E;
    eta_tha(i) = eta_T;
    tmp_array(i) = tmp;
 
    % Compute heat/cold added or removed due to PVT surplus/deficit
    Q_PVT_surplus_actual = max(PVT_output_actual(i) - Total_heat_demand_all_buildings(i), 0); 
    Q_PVT_deficit_actual = max(Total_heat_demand_all_buildings(i) - PVT_output_actual(i), 0); 

    % Compute total heat/cold added or removed from aquifer 1
    Q_added_aquifer_1_actual = thermal_efficiency1*(Q_PVT_surplus_actual*eta_hex - Q_PVT_deficit_actual) + Total_cold_demand_all_buildings(i)*(1-eta_hex)^2;
    % Compute total heat/cold added or removed from aquifer 2
    Q_added_aquifer_2_actual = thermal_efficiency2* (Total_cold_demand_all_buildings(i)) - Q_PVT_deficit_actual*(1-eta_hex)^2;
    
    % Convert Q_added from kWh to Joules
    Q_added_aquifer_1_joules_actual = Q_added_aquifer_1_actual * 3.6e6;
    Q_added_aquifer_2_joules_actual = Q_added_aquifer_2_actual * 3.6e6;

    % Update aquifer temperatures using the Joules values
    delta_T_aquifer_1_actual = Q_added_aquifer_1_joules_actual / (C_w * m_1);
    delta_T_aquifer_2_actual = Q_added_aquifer_2_joules_actual / (C_w * m_2);

    T_aquifer_1_actual(i+1) = T_aquifer_1_actual(i) + delta_T_aquifer_1_actual;
    T_aquifer_2_actual(i+1) = T_aquifer_2_actual(i) + delta_T_aquifer_2_actual;

    % Convert temperatures from Kelvin to Celsius for plotting
    T_aquifer_1_celsius_actual(i) = T_aquifer_1_actual(i) - 273.15;
    T_aquifer_2_celsius_actual(i) = T_aquifer_2_actual(i) - 273.15;
    T_avg_actual = (T_aquifer_1_celsius_actual + T_aquifer_2_celsius_actual)/2;
end

% Calculate and display the maximum temperature difference for each aquifer
maxDiff1_actual = max(T_aquifer_1_celsius_actual) - min(T_aquifer_1_celsius_actual);
maxDiff2_actual = max(T_aquifer_2_celsius_actual) - min(T_aquifer_2_celsius_actual);

% Calculate the difference between T_outflow and T_inflow
T_diff_array = T_outflow_array - T_inflow_array;

%% Calculate CO2 mitigation
% Total heat delivered in kWh
energy_content_ng = 10; % Energy content of natural gas in kWh/m^3
co2_emission_factor_ng = 2.204; % CO2 emissions factor for natural gas in kg CO2 per kWh
combustion_efficiency_ng = 0.9; % Combustion efficiency of natural gas (e.g., 90%)

% Adjust energy content for combustion efficiency
effective_energy_content_ng = energy_content_ng * combustion_efficiency_ng;

% Calculate the equivalent natural gas volume
equivalent_ng_volume = Total_annual_heat_demand / effective_energy_content_ng;

% Calculate CO2 emissions for the equivalent natural gas volume
co2_emissions = equivalent_ng_volume * co2_emission_factor_ng;

% Display the results
fprintf('Equivalent Natural Gas Volume: %f m^3\n', equivalent_ng_volume);
fprintf('CO2 Emissions Avoided: %.0f kg\n', co2_emissions);

%% Efficiency PVT module
% Time loop
for i = 1:(length(hours_year))
    % Update T_i for the current iteration
    if i == 1
        T_inflow = T_initial_2 - 273.15;  
        tmp = 0;
    else
        T_inflow = T_outflow - (T_outflow - T_aquifer_2_celsius_actual(i-1)) * (1-eta_hex);
        T_i = T_inflow; 
        tmp = (T_i - T_aquifer_2_celsius_actual(i-1)) ./ I(i);
    end

    % Second Order Thermal Efficiency Equation
    eta_T = max(min(etao_tha - a1a * tmp - a2a * tmp.^2, 1), 0);
    % Electrical Efficiency Equation
    eta_E = max(min(etao_ele - e1a * tmp, 0.1286), 0);
    % Instantaneous thermal power at normal incidence
    Qa = eta_T .* A_g .* I(i);
    % Instantaneous electrical power at normal incidence
    Qe = eta_E .* A_g .* I(i);

    % Update T_outflow based on T_inflow and Qa
    T_outflow = T_inflow + Qa / (C_w * dm_f);
    PVT_output_actual(i) = Qa * total_PVT_modules / 1000; % In Kw per hour for the entire year
    % Store T_inflow and T_outflow values in the arrays
    T_inflow_array(i) = T_inflow;
    T_outflow_array(i) = T_outflow;
    Qa_array(i) = Qa;
    Qe_array(i) = Qe;
    eta_E_array(i) = eta_E;
    eta_tha(i) = eta_T;
    tmp_array(i) = tmp;
end

% Constants PVT module
eta_elec_STC = 0.18;     % Example: 18% electrical efficiency at STC
beta = 0.0045;           % Example: -0.45% per °C
T_STC = 25;              % Standard Test Condition temperature in °C
NOCT = 48;
G_STC = 1000;

T_in = T_inflow_array; 
Two = T_outflow_array; 

T_pv_normal = Tam + I/G_STC * (NOCT-20);
Tc = (30 + 0.0175 * (I - 300) + 1.14 * (Tam - 25)) + 0.5 * (T_in + Two) - Tam;
    
eta_elec = eta_elec_STC - beta * (T_pv_normal - T_STC);
eta_E = eta_elec_STC - beta * (Tc - T_STC);
eta_E = eta_E';

for i = 1:length(hours_year)
    if I(i) == 0
        eta_elec(i) = 0;  % Set electrical efficiency to zero when irradiance is zero
        eta_E(i) = 0; 
    else
        % Calculate eta_elec and eta_elec1 normally
        deltaT = T_pv_normal(i) - T_STC; % Change in temperature from STC
        eta_elec(i) = eta_elec_STC - beta * deltaT; % Electrical efficiency calculation

        deltaT1 = Tc(i) - T_STC;
        eta_E(i) = eta_elec_STC - beta * deltaT1; 
    end
end                                                       

% Total solar irradiance received by the PVT system in kWh/m^2/year
total_I_sun = sum(mean(WEATHER_output.Irr,2))/1000; 
% Average thermal efficiency
etaSum = sum(eta_tha)/length(I1);
% Average thermal efficiency
etaE = sum(eta_E)/length(I1);

% change if kWh or kWh/m2
% Yearly total thermal output in kWh/m2
total_thermal_output = etaSum*total_I_sun;
% Yearly total electrical output in kWh/m2
total_electrical_output = etaE*total_I_sun;
pvt_collector_output.eta_tha=eta_tha;
pvt_collector_output.etaSum=etaSum;
pvt_collector_output.eta_E=eta_E;
pvt_collector_output.etaE=etaE;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.total_thermal_output=total_thermal_output;
pvt_collector_output.total_electrical_output=total_electrical_output;

%% Plot results
period=WEATHER_output.Period;
Irr = mean(WEATHER_output.Irr,2);
ambient_temperature = WEATHER_output.ambient_temperature;

time = length(Irr);
if time <= 24
    x = 1:time;
    % Irradiance
    f = figure(33); hold on;
 
    yyaxis right
    fill([x, fliplr(x)],[Irr; zeros(time,1)],'y','LineStyle','none')
    alpha(0.33)
    ylabel('Irradiance [W/m^2]');
    f.CurrentAxes.YColor = 'k';
    
    yyaxis left
    plot(x,ambient_temperature,'b-','LineWidth',2);
    f.CurrentAxes.YColor = 'k';
    hold off

    xlabel('Time [Hours]'); ylabel('Temperature [^oC]');
    legend({'T_{amb}','Irradiance'},'Location','NorthEast',...
        'FontSize',14);
    xlim([1 time])
    SetFigureDefaults(f)
else
    xx = 1:time;

% Xticks settings (per hour data)
    months = 24*[31,28,31,30,31,30,31,31,30,31,30,31];  
    ticks = 24*[15.5 31 45 59 74.5 90 105 120 135.5 151 166 181 197 212 227.5 ...
        243 258 273 288 303 317 332 350]-(period(1)+sum(months(1:period(2)-1)));
    month_labelss = {'Jan','','Feb','','Mar','','Apr','','May','','Jun',...
            '','Jul','','Aug','','Sep','','Oct','','Nov','','Dec'};
end

f = figure(1);
% Yearly heat and cold demand for all buildings
plot(xx, normalized_heating_demand_yearly,'DisplayName','Heating Demand');
hold on
plot(xx, normalized_cooling_demand_yearly,'DisplayName','Cooling Demand');
text(2000, 0.9*max(normalized_heating_demand_yearly), sprintf('Total Yearly Heat Demand: %.2f kWh', Total_yearly_heat_demand), 'Color', 'b'); % Display total yearly heat demand
text(2000, 0.8*max(normalized_heating_demand_yearly), sprintf('Total Yearly Cold Demand: %.2f kWh', Total_yearly_cold_demand), 'Color', 'b'); % Display total yearly cold demand
ylabel('Demand (kWh)')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Theoretical yearly Heat and Cold Demand for All Buildings');
xlim([1 time])
legend();
SetFigureDefaults(f)

f = figure(2);
plot(xx, Injected_heat_cooling_season, 'DisplayName', 'Injected Heat (cooling season)');
hold on;
plot(xx, Extracted_cold_cooling_season, 'DisplayName', 'Extracted Cold (cooling season)');
plot(xx, Injected_cold_heating_season, 'DisplayName', 'Injected Cold (heating season)');
plot(xx, Extracted_heat_heating_season, 'DisplayName', 'Extracted Heat (heating season)');
% Add text annotations directly using positions and formatted numbers
y_max = max([max(Injected_heat_cooling_season), max(Extracted_cold_cooling_season), max(Injected_cold_heating_season), max(Extracted_heat_heating_season)]);
text(1, y_max, sprintf('Sum Injected Heat: %.1e kWh', Total_injected_heat), 'Color', 'black');
text(1, y_max - 0.05 * y_max, sprintf('Sum Extracted Cold: %.1e kWh', Total_extracted_cold), 'Color', 'black');
text(1, y_max - 0.1 * y_max, sprintf('Sum Injected Cold: %.1e kWh', Total_injected_cold), 'Color', 'black');
text(1, y_max - 0.15 * y_max, sprintf('Sum Extracted Heat: %.1e kWh', Total_extracted_heat), 'Color', 'black');
ylabel('Energy (kWh)');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Yearly Distribution of estimated ATES Inputs');
xlim([1 time])
legend();
SetFigureDefaults(f)

% Calculate total normalized annual cold demand
% Total_annual_cold_demand = sum(Total_cold_demand_all_buildings);
% Yearly heat demand of all buildings
f = figure(3);
% Yearly heat and cold demand for all buildings
subplot(2,1,1)
hold on
plot(xx, Total_heat_demand_all_buildings,'DisplayName','Heating Demand');
ylabel('Demand (kWh)')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Yearly heat demand of all buildings');
xlim([1 time])
legend();
SetFigureDefaults(f)
subplot(2,1,2)
plot(xx, Total_cold_demand_all_buildings,'DisplayName','Cooling Demand');
ylabel('Demand (kWh)')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Yearly cooling demand of all buildings');
xlim([1 time])
legend();
SetFigureDefaults(f)

f = figure(4);
hold on;
plot(xx, T_aquifer_1_celsius_actual, 'r', 'DisplayName', 'Reservoir 1');
plot(xx, T_aquifer_2_celsius_actual, 'b', 'DisplayName', 'Reservoir 2');
plot(xx, T_avg_actual, 'g', 'DisplayName', 'Average');
ylabel('Temperature (°C)');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Temperature trend of aquifers throughout the year');
xlim([1 time])
legend();
SetFigureDefaults(f)
% Display temperatures at the start of the year
startDay = xx(1);
startTemp1_actual = T_aquifer_1_celsius_actual(1);
startTemp2_actual = T_aquifer_2_celsius_actual(1);
startTempAvg_actual = T_avg_actual(1);
text(startDay, startTemp1_actual, sprintf('%.2f °C', startTemp1_actual), 'Color', 'r', 'VerticalAlignment', 'bottom');
text(startDay, startTemp2_actual, sprintf('%.2f °C', startTemp2_actual), 'Color', 'b', 'VerticalAlignment', 'top');
text(startDay, startTempAvg_actual, sprintf('%.2f °C', startTempAvg_actual), 'Color', 'g', 'VerticalAlignment', 'middle');
% Display temperatures at the end of the year
endDay = xx(end);
endTemp1_actual = T_aquifer_1_celsius_actual(end);
endTemp2_actual = T_aquifer_2_celsius_actual(end);
endTempAvg_actual = T_avg_actual(end);
text(endDay, endTemp1_actual, sprintf('%.2f °C', endTemp1_actual), 'Color', 'r', 'VerticalAlignment', 'bottom');
text(endDay, endTemp2_actual, sprintf('%.2f °C', endTemp2_actual), 'Color', 'b', 'VerticalAlignment', 'top');
text(endDay, endTempAvg_actual, sprintf('%.2f °C', endTempAvg_actual), 'Color', 'g', 'VerticalAlignment', 'middle');

% For Reservoir 1
[minTemp1_actual, minTemp1Index] = min(T_aquifer_1_celsius_actual);
[maxTemp1_actual, maxTemp1Index] = max(T_aquifer_1_celsius_actual);
midPoint1Index = round((minTemp1Index + maxTemp1Index) / 2);

plot(xx(minTemp1Index), minTemp1_actual, 'ro', 'MarkerFaceColor', 'r');
plot(xx(maxTemp1Index), maxTemp1_actual, 'ro', 'MarkerFaceColor', 'r');
line([xx(midPoint1Index), xx(midPoint1Index)], [minTemp1_actual, maxTemp1_actual], 'Color', 'r', 'LineStyle', '-');
line([xx(minTemp1Index), xx(midPoint1Index)], [minTemp1_actual, minTemp1_actual], 'Color', 'r', 'LineStyle', ':');
line([xx(maxTemp1Index), xx(midPoint1Index)], [maxTemp1_actual, maxTemp1_actual], 'Color', 'r', 'LineStyle', ':');

% Display the maximum temperature difference for Reservoir 1
text(xx(midPoint1Index), (minTemp1_actual + maxTemp1_actual) / 2, sprintf('%.2f °C', maxDiff1_actual), 'Color', 'r', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% For Reservoir 2
[minTemp2_actual, minTemp2Index] = min(T_aquifer_2_celsius_actual);
[maxTemp2_actual, maxTemp2Index] = max(T_aquifer_2_celsius_actual);
midPoint2Index = round((minTemp2Index + maxTemp2Index) / 2);

plot(xx(minTemp2Index), minTemp2_actual, 'bo', 'MarkerFaceColor', 'b');
plot(xx(maxTemp2Index), maxTemp2_actual, 'bo', 'MarkerFaceColor', 'b');
line([xx(midPoint2Index), xx(midPoint2Index)], [minTemp2_actual, maxTemp2_actual], 'Color', 'b', 'LineStyle', '-');
line([xx(minTemp2Index), xx(midPoint2Index)], [minTemp2_actual, minTemp2_actual], 'Color', 'b', 'LineStyle', ':');
line([xx(maxTemp2Index), xx(midPoint2Index)], [maxTemp2_actual, maxTemp2_actual], 'Color', 'b', 'LineStyle', ':');

% Display the maximum temperature difference for Reservoir 2
text(xx(midPoint2Index), (minTemp2_actual + maxTemp2_actual) / 2, sprintf('%.2f °C', maxDiff2_actual), 'Color', 'b', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% Plotting and capturing handles for the lines
h1_actual = plot(xx, T_aquifer_1_celsius_actual, 'r', 'DisplayName', 'Aquifer 1');
h2_actual = plot(xx, T_aquifer_2_celsius_actual, 'b', 'DisplayName', 'Aquifer 2');
h_avg_actual = plot(xx, T_avg_actual, 'g', 'DisplayName', 'Average');
legend([h1_actual, h2_actual, h_avg_actual], 'Location', 'northwest');
grid on;
hold off;

% Plotting T_inflow and T_outflow in separate subplots
f = figure(5); 
subplot(2,1,1);
plot(xx, T_inflow_array, 'b-', 'DisplayName', 'T_{inflow}');
ylabel('Temperature (°C)');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('PVT inflow temperature throughout the Year');
xlim([1 time])
legend();
SetFigureDefaults(f)

subplot(2,1,2);
plot(xx, T_outflow_array, 'r-', 'DisplayName', 'T_{outflow}');
ylabel('Temperature (°C)');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('PVT outflow temperature throughout the Year');
xlim([1 time])
legend();
SetFigureDefaults(f)

% Plot the difference
f = figure(7);
plot(xx, T_diff_array, 'DisplayName', 'T_{diff}'); % 'm-' specifies a magenta line for the plot
ylabel('Temperature Difference (°C)');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Difference between outflow and inflow temperatures');
xlim([1 time])
legend()
SetFigureDefaults(f)

% Plotting Qa and Qe
f = figure(8); 
subplot(2,1,1); 
plot(xx, Qa_array,'DisplayName', 'P_{thermal}');
ylabel('Thermal Power (W)');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Thermal Power per PVT module');
xlim([1 time])
legend()
SetFigureDefaults(f)

subplot(2,1,2);
plot(xx, Qe_array, 'r-','DisplayName', 'P_{electrical}');
ylabel('Electrical Power (W)');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Electrical Power per PVT module');
xlim([1 time])
legend()
SetFigureDefaults(f)

f = figure(12);
subplot(2,1,1); 
plot(xx, eta_elec,'DisplayName', '\eta_{electrical, normal}');
ylabel('Electric efficiency');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Electric efficiency per module');
xlim([1 time])
legend()
SetFigureDefaults(f)

subplot(2,1,2);
plot(xx, eta_E, 'r-','DisplayName', '\eta_{electrical}');
ylabel('Electrical efficiency');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
title('Electrical efficiency per module');
xlim([1 time])
legend()
SetFigureDefaults(f)

%%
function SetFigureDefaults(f)
    % Set several figure characteristics for correct representation
    % figure object to be modified

    f.CurrentAxes.FontSize = 14;
    f.CurrentAxes.Box = 'on';
    f.CurrentAxes.XGrid = 'on';
    f.CurrentAxes.YGrid = 'on';
    f.CurrentAxes.XLabel.FontSize = 16;
    f.CurrentAxes.YLabel.FontSize = 16;
    f.CurrentAxes.Title.FontSize = 16;
end

%% Recommendations Based on Simulation Results
% target_avg_temp_actual = (T_initial_1 + T_initial_2);  % Convert from K to °C
avg_temp_start_actual = mean([T_aquifer_1_actual(1), T_aquifer_2_actual(1)]);
avg_temp_end_actual = mean([T_aquifer_1_actual(end), T_aquifer_2_actual(end)]);
ending_avg_temp_actual = mean([T_aquifer_1_actual(end), T_aquifer_2_actual(end)]) - 273.15; % Convert from K to °C
ending_temp_1_actual = T_aquifer_1_actual(end) - T_aquifer_1_actual(1);
ending_temp_2_actual = T_aquifer_2_actual(end) - T_aquifer_2_actual(1);

% Calculate maximum temperature fluctuation in each aquifer
max_temp_fluctuation_1_actual = max(T_aquifer_1_celsius_actual) - min(T_aquifer_1_celsius_actual);
max_temp_fluctuation_2_actual = max(T_aquifer_2_celsius_actual) - min(T_aquifer_2_celsius_actual);

% Calculate minimum temperature difference between aquifers
% min_temp_diff_between_aquifers_actual = min(abs(T_aquifer_1_celsius_actual - T_aquifer_2_celsius_actual));

% Start recommendations
% recommendationMade_actual = false;  % Initialize the recommendation flag at the beginning
% sizeAdjustmentRecommended_actual = false;  % Initialize the size adjustment recommendation flag

% Check if temperatures are within realistic range
if any(T_aquifer_1_celsius_actual > Max_temp) || any(T_aquifer_2_celsius_actual < Min_temp)
    fprintf('Warning: aquifer temperatures are not within the realistic range!\n');
    return;  % Exit the script early
end

% Display the results
fprintf('Ending Average Temperature of the aquifer: %.2f °C\n', ending_avg_temp_actual);
fprintf('Maximum Temperature Fluctuation in aquifer 1: %.2f °C\n', max_temp_fluctuation_1_actual);
fprintf('Maximum Temperature Fluctuation in aquifer 2: %.2f °C\n', max_temp_fluctuation_2_actual);

% Recommendations based on average temperature
if avg_temp_end_actual > avg_temp_start_actual
    if avg_temp_end_actual > avg_temp_start_actual + Average_temp_diff  
        disp('The average temperature has increased and is above the desired level.');
        disp('Consider decreasing the number of PVT modules');
    else
        disp('The average temperature has increased but is still below the desired level.')
        disp('No immediate adjustments for amount of PVT modules necessary.');
    end
elseif avg_temp_end_actual < avg_temp_start_actual
    if avg_temp_end_actual < avg_temp_start_actual - Average_temp_diff
        disp('The average temperature has decreased and is below the desired level.');
        disp('Consider increasing the number of PVT modules');
    else
        disp('The average temperature has decreased but is still above the desired level.');
        disp('No immediate adjustments for amount of PVT modules necessary.');
    end
end

sizeAdjustmentRecommended_actual = false;

% Check for temperature fluctuation in Reservoir 1
if max_temp_fluctuation_1_actual > Max_temp_diff
    fprintf('Aquifer 1 fluctuates too much in temperature. Consider increasing its size.\n');
    sizeAdjustmentRecommended_actual = true;
end

% Check for low temperature fluctuation in Reservoir 1
if max_temp_fluctuation_1_actual < Min_temp_diff2
    fprintf('Aquifer 1 temperature fluctuation is less than the minimum required. Consider decreasing its size.\n');
    sizeAdjustmentRecommended_actual = true;
end

% Check for ending temperature differences in Reservoir 1
if ending_temp_1_actual > Ending_temp_diff_1
    fprintf('Aquifer 1 has warmed up more than is allowed. Consider increasing its size or decreasing the number of PVT modules.\n');
%     sizeAdjustmentRecommended_actual = true;
elseif ending_temp_1_actual < -Ending_temp_diff_1
    fprintf('Aquifer 1 has cooled down more than is allowed. Consider increasing its size.\n');
%     sizeAdjustmentRecommended_actual = true;
elseif ~sizeAdjustmentRecommended_actual
    fprintf('Aquifer 1 is within the expected temperature range.\n');
end

sizeAdjustmentRecommended_actual = false;

% Check for temperature fluctuation in Reservoir 2
if max_temp_fluctuation_2_actual > Max_temp_diff
    fprintf('Aquifer 2 fluctuates too much in temperature. Consider increasing its size.\n');
    sizeAdjustmentRecommended_actual = true;
end

% Check for low temperature fluctuation in Reservoir 2
if max_temp_fluctuation_2_actual < Min_temp_diff2
    fprintf('Aquifer 2 temperature fluctuation is less than the minimum required. Consider decreasing its size.\n');
    sizeAdjustmentRecommended_actual = true;
end

% Check for ending temperature differences in Reservoir 2
if ending_temp_2_actual > Ending_temp_diff_2
    fprintf('Aquifer 2 has warmed up more than expected. Consider increasing its size.\n');
%     sizeAdjustmentRecommended_actual = true;
elseif ending_temp_2_actual < -Ending_temp_diff_2
    fprintf('Aquifer 2 has cooled down more than expected. Consider increasing its size.\n');
%     sizeAdjustmentRecommended_actual = true;
elseif ~sizeAdjustmentRecommended_actual
    fprintf('Aquifer 2 is within the expected temperature range.\n');
end                                          
end