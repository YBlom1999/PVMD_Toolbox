function [electricityconsumption,boilerW] = heatpumpint(TOOLBOX_input,Size_HP)
% heatpumpint calculates the electricity consumption of the residential
% building
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
% Size_HP: double
%   
%
% Returns
% -------
% electricityconsumption : double
%   The electricity consumption over the whole year
% boilerW: double
%   The energy produced by the boiler.
%
% Developed by Tekin Kodzhabash, Implemented by: Y. Blom
%Loading the input
[weather_data, ~] = load_meteonorm_data(TOOLBOX_input);
T_amb = weather_data(:, 9); 
%data from npro
[~,~,data_folder] = get_folder_structure;
mainexcel = fullfile(data_folder, 'Financial','NproData',TOOLBOX_input.FinancialPart.NproFile);
heatsheet = 'Heat';
spaceheating_npro = xlsread(mainexcel, heatsheet, 'C18:C8777');
dailydhw= xlsread(mainexcel,heatsheet,'H18:H382');
coldsheet='Cold';
spacecooling_npro=xlsread(mainexcel, coldsheet, 'C18:C8777');

boiler=zeros(8760,1);
spaceheating_hp=zeros(8760,1);
dhw_hp=zeros(8760,1);

%fixing space heating demand
index = (T_amb > 15); % Create a logical index where T_amb_ams is more than 15
spaceheating_demand= spaceheating_npro; % Initialize 'powersupply_ASHP_radiator_ams_fixed' as a copy of 'powersupply_ASHP_radiator_ams'
spaceheating_demand(index) = 0; % Set the values to zero for the positions where T_amb_ams is more than 15

max_spaceheating=max(spaceheating_demand);

%fixing space cooling demand
indexc=(T_amb<25); %Create a logical index where T_amb_ams is less than 25
spacecooling_demand=spacecooling_npro;
spacecooling_demand(indexc)=0; %Set the values to zero for the positions where T_amb_ams is less than 25

max_spacecooling=max(spacecooling_demand);

%creating fixed daily dhw
days_in_year = 365;
hours_in_day = 24;
dhw_fixed_annual = zeros(days_in_year * hours_in_day, 1);% Initialize dhw_fixed_annual as a column vector with zeros
hour_indices = (11:15);% Define the indices for hours 11 to 15
for day = 1:days_in_year % Loop through each day in a year and set the values in dhw_fixed_annual
    % Calculate the starting and ending indices for the current day
    start_index = (day - 1) * hours_in_day + hour_indices(1);
    end_index = (day - 1) * hours_in_day + hour_indices(end);

    % Set the values for hours 11 to 14 for the current day
    dhw_fixed_annual(start_index:end_index) = dailydhw(day) / numel(hour_indices);
end

%To consider the boiler action when the total energy to be supplied is
%higher than the heat pump capacity
for i=1:8760
    if spaceheating_demand(i)+dhw_fixed_annual(i)>Size_HP %heat pump capacity is less than what needs to be supllied for heating
        if spaceheating_demand(i)>Size_HP %heat pump capacity is less than space heating demand
            spaceheating_hp(i)=Size_HP; %heat pump gives everything to space heating and disregard dhw
            boiler(i)=dhw_fixed_annual(i)+spaceheating_demand(i)-spaceheating_hp(i) ;%amount of energy boiler needs to give
            dhw_hp(i)=0; %no energy for dhw because it is supplied by boiler
        else %space heating can be supplied by heat pump completeley and some of dhw can be supplied by heat pump
            spaceheating_hp(i)=spaceheating_demand(i);
            boiler(i)=dhw_fixed_annual(i)+spaceheating_hp(i)-Size_HP;
            dhw_hp(i)=Size_HP-spaceheating_hp(i);
        end
    else %everything can be supplied by the heat pump
        spaceheating_hp(i)=spaceheating_demand(i);
        dhw_hp(i)=dhw_fixed_annual(i);
    end   
end

% To eliminate cooling and the heating at the same time, space cooling and
% space heating cannot happen at the same time anyway.
for i=1:8760
    if spacecooling_demand(i)>0
        boiler(i)=dhw_fixed_annual(i);
        dhw_hp(i)=0;
    end
end

% To eliminate the cooling more than the heat pump capacity
spacecooling_demand_real=spacecooling_demand;
for i=1:8760
    if spacecooling_demand(i)>Size_HP
        spacecooling_demand(i)=Size_HP;
    end

end


%COP calculations with radiator supply for spaceheating
COP_ASHP_radiator= 0.85*cashp(radiator(T_amb),T_amb);% 0.85 is the correction factor
%COP calculation for dhw 50 degrees is the sink side temperature of the
%water supply
COP_ASHP_dhw=0.85*cashp(50,T_amb);

%EER calculation for cooling
EER_ref=5.8115+0.0014.*(T_amb-10).^2-0.154.*(T_amb-10);

%power calculations
powersupply_ASHP_radiator=spaceheating_hp./COP_ASHP_radiator;
powersupply_ASHP_dhw=dhw_hp./COP_ASHP_dhw;
powersupply_ASHP_cooling=spacecooling_demand./EER_ref;

%electricity load
electricitysheet = 'Electricity';
plugload = xlsread(mainexcel, electricitysheet, 'C18:C8777'); %kW


%total electricity consumption

electricityconsumption=(plugload+powersupply_ASHP_dhw+ powersupply_ASHP_radiator+powersupply_ASHP_cooling)*1000; %W
boilerW=boiler*1000; %changin boiler to kwH from Wh
end



function [radiator] = radiator(T)
% radiator calculates the temperature difference for a radiator
%
% Parameters
% ----------
% T : double
%   Ambient temperature
%   
%
% Returns
% -------
% radiator : double
%   Temperature difference
%
% Developed by Tekin Kodzhabash, Implemented by: Y. Blom
radiator = 40-T;
end

function [COP_ASHP] = cashp(T_sink,T_source)
% cashp calculates the COP of the ASHP
%
% Parameters
% ----------
% T_sink : double
%   Temperature of the sink
% T_source : double
%   Temperature of the source
%   
%
% Returns
% -------
% COP_ASHP: double
%   The coefficient of performance
%
% Developed by Tekin Kodzhabash, Implemented by: Y. Blom
%COP of ASHP
COP_ASHP = 6.08-0.09*(T_sink-T_source)+0.0005*(T_sink-T_source).^2;
end