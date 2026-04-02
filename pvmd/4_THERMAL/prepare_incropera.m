function [pvModule,weather] = prepare_incropera(TOOLBOX_input, WEATHER, CELL_output,cell_i)


pvModule.tilt = TOOLBOX_input.Scene.module_mounting.ModTilt;
pvModule.azimuth = TOOLBOX_input.Scene.module_mounting.ModAzimuth;
pvModule.albedo = TOOLBOX_input.Scene.module_mounting.Albedo;

CR = TOOLBOX_input.Scene.module_mounting.CellRows;
CC = TOOLBOX_input.Scene.module_mounting.CellColumns;
CS = TOOLBOX_input.Scene.module_mounting.CellSpacing*1e-2;
ES = TOOLBOX_input.Scene.module_mounting.EdgeSpacing*1e-2;
CL = TOOLBOX_input.Scene.module_mounting.CellLength*1e-2;
CW = TOOLBOX_input.Scene.module_mounting.CellWidth*1e-2;
pvModule.length = CR*CL + (CR-1)*CS + 2*ES;
pvModule.width = CC*CW + (CC-1)*CS + 2*ES;
pvModule.layers = TOOLBOX_input.thermal.layers;


weather_data = load_climate_data(TOOLBOX_input.irradiation);
weather.dhi = weather_data(:,8);
weather.dni = weather_data(:,7);
weather.ghi = weather_data(:,11);
weather.tAmb = weather_data(:,9);
weather.windSpeed = weather_data(:,10);
weather.solarAzi = weather_data(:,5);
weather.solarAlti = weather_data(:,6);
weather.tGnd = (weather.tAmb+273.15)+(0.015*(1-pvModule.albedo).*weather.ghi-0.7).*exp(-0.09*weather.windSpeed)-273.15;%(Celsius)
weather.tSky = 0.0552*(weather.tAmb+273.15).^1.5-273.15;%(Celsius) Obtained from the Solar Energy Book
qGen = zeros(TOOLBOX_input.thermal.Nlayers,length(WEATHER.Irr));
assignment = TOOLBOX_input.thermal.assignment_layers;
for i = 2:length(assignment)-1
    if sum(i == CELL_output.CELL_FRONT.Absmat_ind)
        qGen(assignment(i),:) = qGen(assignment(i),:)+WEATHER.A(:,cell_i,i-1)'*(1-TOOLBOX_input.thermal.Efficiency/100);
    else
        qGen(assignment(i),:) = qGen(assignment(i),:)+WEATHER.A(:,cell_i,i-1)';
    end
end
weather.qGen = qGen;
end

function weather_data = load_climate_data(Irradtion)
%Read the meteonorm file
%The code is based on the load_meteonorm_data
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

%---- Load weather file
[~,~,data_folder] = get_folder_structure;
locations_dir = fullfile(data_folder,'Weather','Locations');
climate_file = fullfile(locations_dir,Irradtion.climateFile);
load(climate_file,'weather_data');

%---- Filter weather data based on period selected by user
% The year should be selected not set randomly to 2021
init_day = Irradtion.init_day;
init_month = Irradtion.init_month;
end_day = Irradtion.end_day;
end_month = Irradtion.end_month;
year = Irradtion.year_choice;
[init_day,init_month,end_day,end_month] = ...
    clean_period(init_day,init_month,end_day,end_month,year);
init_cell = (sum(eomday(year,1:init_month-1))+ init_day-1)*24+1;
end_cell = (sum(eomday(year,1:end_month-1))+ end_day-1)*24 + 24;
weather_data = weather_data(init_cell:end_cell,:);

end
