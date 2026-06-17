function [pvt_collector_output] = HP_new_ST_air(~,~, WEATHER_output)
%This function calculates the HP components temperatures, electric power and COP of HP 
% It is integrated with solar_thermal_air
% Ask for one of the coolants & heating_load
% CoolProp open source property generator (C++ library ) => Calculates thermophysical properties for a wide range of fluids, including refrigerants
% Allows to obtain properties such as temperature, pressure, density, enthalpy, entropy, and more for a given fluid at a specific state.
% 02-24
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% MODULE_output : struct
%   Output of this module block
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% pvt_collector_output : struct
%   Simulation results of the PVT module
% TOOLBOX_input : struct
%   Update of TOOLBOX_input by adding simulation parameters of thermal module
% 
% Developed by ZUA.

pyenv;                          % Calling python 3.7 from matlab - use cmd to install coolprop using: 'pip install CoolProp' ---- changed (pyversion)

A_g = 2.6;                      % Gross Collector Area (m2)
dm_f = 0.0844;                  % Fluid Mass Flowrate (kg/m2s)

% Second Order Thermal Efficiency Equation unique coefficients
a1a = 6.48;                     % W/(m2K)                                          
a2a = 0.028;                    % W/(m2K2) 
etao_tha = 0.40;
C_f = 1007;                     % specific heat (J/kg*C) 

I = mean(WEATHER_output.Irr{1},2)';
Tam = [WEATHER_output.ambient_temperature]';

idx = find(I~=0);                % find the indices of non-zero elements in I
if length(idx) > length(I)       % if there are non-zero elements in I
    idx = idx(1:length(I));      % select the non-zero elements
end
I1 = I(idx);                     % select the corresponding elements from I

T_i = Tam + 0.5;                 % Temperature of the fluid entering the collector - assumed (for air-based) (*C)

% preallocate arrays
eta_E = zeros(length(Tam), 1);
Tc = zeros(length(Tam), 1);

% calculate the value that is repeatedly used in the equation for eta_tha
tmp = (T_i - Tam) ./ I;
% Second Order Thermal Efficiency Equation
eta_tha = etao_tha - a1a * tmp - a2a * I .* tmp.^2;
% avoid negative values
eta_tha = max(eta_tha, 0);
% avoid values greater than one
eta_tha = min(eta_tha, 1);
% Instantaneous power at normal incidence
Qa = eta_tha .* A_g .* I;
% fluid outlet temperature
Two = T_i + (Qa ./ (dm_f * C_f));
Tmean = (T_i + Two) ./ 2;
% Total solar irradiance received by the PVT system in kWh/m^2/year
total_I_sun = sum(mean(WEATHER_output.Irr{1},2))/1000; 
% Average thermal efficiency
etaSum = sum(eta_tha)/length(I1);

%% change if kWh or kWh/m2
% Yearly total thermal output in kWh/m2
total_thermal_output = etaSum*total_I_sun*A_g/A_g;
etaE = 0;                               % because it is a solar thermal collector only - just for toolbox
total_electrical_output = 0;
pvt_collector_output.eta_tha=eta_tha;
pvt_collector_output.etaSum=etaSum;
pvt_collector_output.eta_E=eta_E;
pvt_collector_output.etaE=etaE;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.Tmean=Tmean;
pvt_collector_output.total_thermal_output=total_thermal_output;
pvt_collector_output.total_electrical_output=total_electrical_output;

%% Heat pump
%% Prompt the user to select a coolant
% Display coolants
disp('Coolant Options:');
disp('1. R134a');
disp('2. R410a');
disp('3. R407c');
% R134a (Hydrofluorocarbon refrigerant): Widely used in various refrigeration and air conditioning applications.
% R410A (Mixture of Difluoromethane and Pentafluoroethane): Widely used in air-source heat pumps for both heating and cooling applications.
% R407C (Mixture of Difluoromethane, Pentafluoroethane, and Tetrafluoroethane): Used in air-source heat pumps as a replacement for R22.
% User to select a coolant
coolant_prompt = ('Select the number for your preferred coolant (1. R134a, 2. R410a, 3. R407c):');
coolant_dlgtitle = 'Coolant Selection';
% Get user input
coolant_choice = inputdlg(coolant_prompt, coolant_dlgtitle);
% Check if the user entered a value or canceled
if isempty(coolant_choice)
    disp('Coolant selection canceled.');
else
    % Convert input string to a numeric value for coolant
    coolant_number = str2double(coolant_choice);
    % Check if the input is a valid number (1, 2, or 3) for coolant
    if isnan(coolant_number) || ~ismember(coolant_number, [1, 2, 3])
        disp('Invalid input for coolant. Please select 1, 2, or 3.');
    else
        % Update the fluid variable based on the selected coolant
        coolants = {'R134a', 'R410a', 'R407c'};
        fluid = coolants{coolant_number};
        % Display the selected coolant
        disp(['Selected coolant: ', fluid]);
        % User to enter the heat load
        heat_load_prompt = 'Enter the heat load (1 to 10000):';
        heat_load_dlgtitle = 'Heat Load Input';
        heat_load_input = inputdlg(heat_load_prompt, heat_load_dlgtitle);
        % Check if the user entered a value or canceled
        if isempty(heat_load_input)
            disp('Heat load input canceled.');
        else
            % Convert input string to a numeric value for heat load
            Q_load = str2double(heat_load_input);
            % Check if the input is a valid number (between 1 and 10000) for heat load
            if isnan(Q_load) || Q_load < 1 || Q_load > 10000
                disp('Invalid input for heat load. Please enter a value between 1 and 10000.');
            else
                % Use the updated fluid and heat_load variables
                disp(['Perform calculations with ', fluid, ' and heat load ', num2str(Q_load)]);
            end
        end
    end
end

% These can be changed 
T_des = 40;                                     % Desired temperature of room
eta_comp = 0.75;                                % Compressor isentropic efficiency (use any value between 0 and 1)
T_evap_appr = 5;                                % Evaporator approach temperature (temperature between the working fluid leaving the condenser and the ambient temperature)
T_cond_appr = 5;                                % Condenser approach temperature (temperature between the room temperature and the working fluid leaving the condenser)

temp = size(Tam);                               % Determine size of the Tamb array
n = temp(1,2);                                  % Extracting the size of Tamb array
temp = 0;                                       % Defining to be used in evaluating average COP
E_input = 0;                                    % Defining to store cummulative energy input in kWh
T1 = zeros(1,n);                                % Defining T1 and T4 array
T2 = zeros(1,n);                                
T3 = zeros(1,n);                                
P_hp = zeros(1,n);                              % Compressor work
P_hp_new = zeros(1,n); 
Q_room = zeros(1,n);                            % Heat transferred in the condenser to the room heat load
Q_amb = zeros(1,n);                             % Heat transferred in the evaporator to the ambient
COP_hp = zeros(1,n);

% Solve heat pump for each ambient temperature
for i=1:n
    T1(1,i) = Two(i) - T_evap_appr;                                                    % Evaporator temperature in °C
    % T4 = T1 because in the two phase region and at same pressure
    % P1 = py.CoolProp.CoolProp.PropsSI("P","T",T1(1,i) + 273.15,"Q",1,fluid);         % Evaporator pressure in Pa
    T3(1,i) = T_des + T_cond_appr;                                                     % Condenser temperature in °C
    P2 = py.CoolProp.CoolProp.PropsSI("P","T",T3(1,i) + 273.15,"Q",0,fluid);           % Condenser pressure in Pa
    % P3 = P2 as there is no pressure drop
    h1 = py.CoolProp.CoolProp.PropsSI("H","T",T1(1,i) + 273.15,"Q",1,fluid)/1000;      % Enthalpy at point 1 in kJ/kg
    s1 = py.CoolProp.CoolProp.PropsSI("S","T",T1(1,i) + 273.15,"Q",1,fluid)/1000;      % Entropy at point 1 in kJ/kg °C
    s2s = s1;                                                                          % Entropy at point 2 considering isentropic compressor
    h2s = py.CoolProp.CoolProp.PropsSI("H","P",P2,"S",s2s*1000,fluid)/1000;            % Enthalpy at point 2 in kJ/kg considering isentropic compressor
    % h2s is when compressor isentropic efficiency is 100% . In this case compression process is ideal when its 75% then the process is not ideal, and enthalpy after compression is h2a
    h2a = h1 + (h2s - h1)/eta_comp;                                                    % Enthalpy at point 2 in kJ/kg considering compressor isentropic efficiency
    T2(i) = py.CoolProp.CoolProp.PropsSI("T","H",h2a*1000,"P",P2,fluid) - 273.15;      % Temperature after compressions (compressed gas maybe in superheated form after compression)
    P_hp(1,i) = (h2a-h1);                                                              % Compressor work in kW
    h3 = py.CoolProp.CoolProp.PropsSI("H","T",T3(1,i) + 273.15,"Q",0,fluid)/1000;      % Enthalpy at point 3 in kJ/kg
    h4 = h3;                                                                           % Enthalpy at point 4 in kJ/kg
    Q_amb(1,i) = (h1 - h4);                                                            % Heat duty of the heat pump in kW
    eta_pvt_hp = 0.60;
    mf  = dm_f*(Two(1,i) - T_i(1,i))*C_f*eta_pvt_hp/Q_amb(1,i);
    % Q_evap = Q_amb
    Q_room(1,i) = (h2a - h3);                                                          % Heat transfer in the condenser to the room
    % Q_cond = Q_room
    P_hp_new(1,i) = mf*P_hp(1,i);
    COP_hp(1,i) = (Q_room(1,i))/P_hp(1,i);                                             % COP   
    E_input = E_input + (Q_load/(365*24))*P_hp_new(1,i)/Q_amb(1,i);                        % Cummulatively adding the energy input
    temp = temp + COP_hp(1,i);                                                         % Cummulatively adding COP to be used later for average COP calculations
end

pvt_collector_output.P_hp_total(1,i)=P_hp_new(1,i);
pvt_collector_output.COP_hp=COP_hp;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.total_thermal_output=total_thermal_output;
pvt_collector_output.total_electrical_output=total_electrical_output;

COP_avg = temp/n;                                                                      % Average COP
fprintf('Power input for given period is %.4f kWh',E_input);                           % Print daily energy input in kWh
fprintf('\nAveragen COP is %.4f\n',COP_avg);                                           % Print average COP

% Plot
%% Plot results
period=WEATHER_output.Period;
Irr = mean(WEATHER_output.Irr{1},2);
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

 f = figure(10231);
 subplot(1,2,1)
plot(xx,T1);
hold on
plot(xx,T2);
hold on
plot(xx,T3);
plot(xx,Tam);
legend('T_{1,ev-cp} & T_{4,ev-ex}','T_{2,cp-cd}','T_{3,cd-ex}', 'T_{amb}')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Ambient Temperature and HP Component Temperature');
    xlim([1 time])
    ylabel('Temperature [°C]')
        SetFigureDefaults(f)
subplot(1,2,2)
yyaxis left;
plot(xx,P_hp_new)
ylabel('Power input - Compressor work [W]');
yyaxis right;
plot(xx,COP_hp)
ylabel('COP')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('COP & Power input');
    xlim([1 time])
    legend('Compresor work','COP')
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

end