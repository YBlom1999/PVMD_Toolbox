function [weather_data, number_hours] = load_meteonorm_data(TOOLBOX_input)
%LOAD_METEONORM_DATA Read Meteonorm data
% Read a .mat file with meteonorm data for one year for the desired
% location and select the range of interest considering the user input
%
% The columns of the weather_data variable are:
%   column 1 - year
%   column 2 - month
%   column 3 - day
%   column 4 - hour
%   column 5 - sun azimuth (S = 0, W = 90, N = -180, E = -90)
%   column 6 - sun altitude
%   column 7 - DNI
%   column 8 - DHI
%   column 9 - ambient temperature
%   column 10 - wind speed
%   column 11 - GHI
%   column 12 - Relative humidity 
%   column 13 - UVa radiation
%   column 14 - UVb radiation
% The last three columns are only needed for degradation purposes
%
% Use PVMD.muf file as meteonorm template in location: 
% Meteotest\Meteonorm7\CustomFormats\
%
% Note: Meteonorm only works for Windows OS
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
%
% Returns
% -------
% weather_data : double
%   weather data for the selected period%
% number_hours : double
%   length of the weather_data variable
%
% Developed by unknown (E. Garcia?). Commented by A. Alcaniz

%---- Load weather file
[~,~,data_folder] = get_folder_structure;
locations_dir = fullfile(data_folder,'Weather','Locations');
climate_file = fullfile(locations_dir,TOOLBOX_input.irradiation.climateFile);
load(climate_file,'weather_data');

%---- Filter weather data based on period selected by user
% The year should be selected not set randomly to 2021
init_day = TOOLBOX_input.irradiation.init_day;
init_month = TOOLBOX_input.irradiation.init_month;
end_day = TOOLBOX_input.irradiation.end_day;
end_month = TOOLBOX_input.irradiation.end_month;
year = TOOLBOX_input.irradiation.year_choice;
[init_day,init_month,end_day,end_month] = ...
    clean_period(init_day,init_month,end_day,end_month,year);
init_cell = (sum(eomday(year,1:init_month-1))+ init_day-1)*24+1;
end_cell = (sum(eomday(year,1:end_month-1))+ end_day-1)*24 + 24;
weather_data = weather_data(init_cell:end_cell,:);

%---- Display the selected period to the user
init_month_char = month(datetime(2000,init_month,1),'name');
end_month_char = month(datetime(2000,end_month,1),'name');
if init_month == end_month && init_day == end_day
    fprintf('\tSimulating for %s, %d.\n',...
        init_month_char{1},init_day);
else
   fprintf('\tSimulating from %s, %d to %s, %d.\n',...
       init_month_char{1},init_day,...
       end_month_char{1},end_day);
end

%---- Compute the number of hours
number_hours = size(weather_data,1);

end
