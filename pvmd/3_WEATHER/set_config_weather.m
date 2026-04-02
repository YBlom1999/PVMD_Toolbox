function TOOLBOX_input = set_config_weather(TOOLBOX_input)
%SET_CONFIG_WEATHER Get user input if parameters are missing
% Used only for the GUI. Some parameters are first asked to the user, and
% later they are given the proper format for the TOOLBOX_input
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
%
% Returns
% -------
% TOOLBOX_input : struct
%   Simulation parameters

climate_file = select_file('Locations','Choose a climate');
if isempty(climate_file), return, end

[init_day,init_month,end_day,end_month,year] = specify_period();
if isempty(init_day), return, end

[include_horizon,recalculate_horizon,lat,lon,radius,skyline_file] ...
    = choose_horizon_details();
if isempty(include_horizon), return, end

spectra_choice = ChooseSpectrum();
if isempty(spectra_choice), return, end

TOOLBOX_input.runWeatherPart = true;
TOOLBOX_input.irradiation.plot_weather = true;
TOOLBOX_input.irradiation.climateFile = climate_file;
TOOLBOX_input.irradiation.init_day = init_day;
TOOLBOX_input.irradiation.init_month = init_month;
TOOLBOX_input.irradiation.end_day = end_day;
TOOLBOX_input.irradiation.end_month = end_month;
TOOLBOX_input.irradiation.year_choice = year;
TOOLBOX_input.irradiation.spectra_choice = spectra_choice;

TOOLBOX_input.irradiation.include_horizon_reconstruction = include_horizon;
if include_horizon
    TOOLBOX_input.irradiation.plot_skyline = true;
    TOOLBOX_input.irradiation.recalculate_horizon = recalculate_horizon;
    TOOLBOX_input.irradiation.skyline_file = skyline_file;
    if recalculate_horizon
        TOOLBOX_input.irradiation.latitude = lat;
        TOOLBOX_input.irradiation.longitude = lon;
        TOOLBOX_input.irradiation.radius = radius;
    end
end

end


function file_name = select_file(folder_name,prompt)
%SELECT_CLIMATE_FILE Choose one of climate files from the locations folder
%
% Returns
% -------
% climate_file : char
%   Name of the selected climate file

[~,~,data_folder] = get_folder_structure;
files_dir = fullfile(data_folder,'Weather',folder_name);

list = dir([files_dir,filesep,'*.mat']);
A = listdlg('PromptString',prompt,'SelectionMode','single',...
    'ListString',{list.name},'ListSize',[200,100]);
if isempty(A), return, end

file_name = list(A).name;
end


function [init_day,init_month,end_day,end_month,year] = specify_period()
%SPECIFY_PERIOD Ask the user to specify the period of simulation
%
% Returns
% -------
% init_day : double
%   Initial day of the initial month for the simulations
% init_month : double
%   Initial month for the simulations
% end_day : double
%   Final day of the final month for the simulations
% end_month : double
%   Final month for the simulations
% year : double
%   Year of choice for the simulations

prompt = {'Start: Month','Start: Day','End: Month','End: Day','Year of choice'};
default = {'1','1','1','31','2021'};
answer = str2double(inputdlg(prompt,'Specify the period of simulation',...
    1,default));

if isempty(answer), return, end

% Clean up input and assign
% Year is set to 2021 randomly, but it should ask the user which year we
% are in, or at least whether it's a leap year or not
[init_day,init_month,end_day,end_month] = ...
    clean_period(answer(2),answer(1),answer(4),answer(3),answer(5));
year = answer(5);

end


function [include_horizon,recalculate_horizon,lat,lon,radius,...
    skyline_file] = choose_horizon_details()
%CHOOSE_HORIZON_DETAILS Ask the user if the horizon is to be included and
%the needed parameters for it
%
% Returns
% -------
% include_horizon : logical
%   Whether the horizon is to be included or not
% recalculate_horizon : logical
%   Whether the horizon is to be reconstructed or not
% lat : double
%   Latitude for the reconstruction of the horizon
% lon : double
%   Longitude for the reconstruction of the horizon
% radius : double
%   Radius [m] which determines the area selected for the reconstruction
% city : char
%   City where the horizon is reconstructed. For easy identification
% skyline_file : char
%   Selected reconstructed horizon file

question1 = 'Do you want to include the horizon?';
include_horizon = ask_user_yes_no_question(question1);
if include_horizon
    question2 = 'Do you want to calculate the horizon? (No means use an existant one)';
    recalculate_horizon = ask_user_yes_no_question(question2);
    
    if recalculate_horizon
        [lat,lon,radius,city] = get_horizon_reconstruction_params();
        skyline_file = strcat('reconstructed_horizon_',city{1},'.mat');
    else
        skyline_file = select_file('Reconstructed Horizons',...
            'Choose a skyline');
        lat = nan; lon = nan; radius = nan;
    end
else
    recalculate_horizon = nan; skyline_file = nan;
    lat = nan; lon = nan; radius = nan;
end
end


function answer = ask_user_yes_no_question(ques)
% Ask the user a yes/no question through the command window
%
% Parameters
% ----------
% ques : char
%   Question to be asked
%
% Returns
% -------
% answer : logical
%   User answer

Q1 = menu(ques, 'Yes', 'No');
answer = Q1 == 1;
end


function [latitude,longitude,radius,city] = get_horizon_reconstruction_params()
%GET_HORIZON_RECONSTRUCTION_PARAMS Get data needed to construct the horizon
%
% Returns
% -------
% latitude : double
%   Latitude for the reconstruction of the horizon
% longitude : double
%   Longitude for the reconstruction of the horizon
% radius : double
%   Radius [m] which determines the area selected for the reconstruction
% city : char
%   City where the horizon is reconstructed. For easy identification


prompt = {'Latitude','Longitude','Radius [m]','City'};
default = {'52.011753','4.359364','200','Delft'};
answer = inputdlg(prompt,'Specify the horizon reconstruction parameters',...
    1,default);
latitude = str2double(answer(1));
longitude = str2double(answer(2));
radius = str2double(answer(3));
city = answer(4);
end

function spectra_choice = ChooseSpectrum()
%ChooseSpectrum Ask the user to specify the spectrum used for the
%simulation (SMARTS or SBDarts)
%
% Returns
% -------
% spectra_choice : double
%   Choice for the spectra (1 is SMARTS, 2 is SBDarts)

Q ={'Select spectra'};    %ask user
spectra_choice = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'SMARTS (clear sky)','SBDarts (include clouds)'}); %user choice

end