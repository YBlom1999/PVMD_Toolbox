function [pvt_collector_output] = Solar_PV_thermal_different_absorbers(~,~, WEATHER_output)
% Calculates the thermal and PV yield of PVT collector with different absorbers (unique coefficients)
% Thermo photovoltaic collector - pap-I - 5.1 m2
%
% Parameters
% ----------
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% pvt_collector_output : struct
%   Simulation results of the PVT module
% 
% Developed by ZUA.

A_g = 0.64;                     % Gross Collector Area (m2)
dm_f = 0.01;                    % Fluid Mass Flowrate (kg/s)
C_f = 4186;                     % specific heat (J/kg°C) 

% Prompt user to select an absorber
absorber_prompt = 'Select the number for your preferred absorber (1. Direct flow, 2. Oscillatory flow, 3. Serpentine flow, 4. Web flow, 5. Spiral flow, 6. Parallel Serpentine flow, 7. Modified Parallel Serpentine flow):';
absorber_dlgtitle = 'Absorber Selection';
% Get user input
absorber_choice = inputdlg(absorber_prompt, absorber_dlgtitle);

% Check if the user entered a value or canceled
if isempty(absorber_choice)
    disp('Absorber selection canceled.');
else
    % Convert input string to a numeric value for absorber
    absorber_number = str2double(absorber_choice);
    % Check if the input is a valid number for absorber
    if isnan(absorber_number) || ~ismember(absorber_number, [1, 2, 3, 4, 5, 6, 7])
        disp('Invalid input for absorber. Please select a number between 1 and 7.');
    else
        % Update the absorber variables based on the selected absorber
        absorbers = {'Direct flow', 'Oscillatory flow', 'Serpentine flow', 'Web flow', 'Spiral flow', 'Parallel Serpentine flow', 'Modified Parallel Serpentine flow'};
        a1_values = [2.675, 2.31, 1.80, 2.29, 2.619, 2.52, 2.51];                         % Corresponding values for a1a W/(m2K)
        etao_tha_values = [0.4641, 0.4124, 0.3235, 0.4215, 0.5012, 0.4812, 0.4833];       % Corresponding values for etao_tha
        e1a_values = [0.0211, 0.0269, 0.0615, 0.025, 0.00769, 0.01538, 0.01384];          % Corresponding values for e1a W/(m2K)
        etao_ele_values = [0.1185, 0.1181, 0.1157, 0.1183, 0.1194, 0.1189, 0.1190];       % Corresponding values for etao_ele

        absorber = absorbers{absorber_number};
        a1a = a1_values(absorber_number);
        etao_tha = etao_tha_values(absorber_number);
        e1a = e1a_values(absorber_number);
        etao_ele = etao_ele_values(absorber_number);

        % Display the selected absorber and corresponding values
        disp(['Selected absorber: ', absorber]);        
    end
end

I = mean(WEATHER_output.Irr,2)';
Tam = [WEATHER_output.ambient_temperature]';

idx = find(I~=0);                % find the indices of non-zero elements in I
if length(idx) > length(I)       % if there are non-zero elements in I
    idx = idx(1:length(I));      % select the non-zero elements
end
I1 = I(idx);                     % select the corresponding elements from G

T_i = Tam + 2.5;                 % Temperature of the fluid entering the collector - assumed (*C)

% calculate the value that is repeatedly used in the equation for eta_tha
tmp = (T_i - Tam) ./ I;
% Second Order Thermal Efficiency Equation
eta_tha = etao_tha - a1a * tmp;
% avoid negative values
eta_tha = max(eta_tha, 0);
% avoid values greater than one
eta_tha = min(eta_tha, 1);
% Electrical Efficiency Equation
eta_E = etao_ele - e1a * tmp;
% avoid negative values
eta_E = max(eta_E, 0);
% avoid values greater than one
eta_E = min(eta_E, 0.15);
% Instantaneous thermal power at normal incidence
Qa = eta_tha .* A_g .* I;
% Instantaneous electrical power at normal incidence
Qe = eta_E .* A_g .* I;
% fluid outlet temperature
Two = T_i + (Qa ./ (dm_f * C_f));
Tmean = (T_i + Two) ./ 2;
% Cell temperature
NOCT = 48;                       % typical module at 48°C (best module operated at a NOCT of 33°C, the worst at 58°C)
S = I*0.1;                       % insolation in mW/cm^2
Tc = Tam + S*(NOCT - 20)/80;
% Tc_pv = 30 + 0.0175*(I - 300) + 1.14*(Tam - 25);
% Tc = Tc_pv + ((T_i + Two)/2 - Tam);

% Total solar irradiance received by the PVT system in kWh/m^2/year
total_I_sun = sum(mean(WEATHER_output.Irr,2))/1000; 
% Average thermal efficiency
etaSum = sum(eta_tha)/length(I1);
% Average electrical efficiency
etaE = sum(eta_E)/length(I1);

%% change if kWh or kWh/m2
% Yearly total thermal output in kWh/m2
total_thermal_output = etaSum*total_I_sun*A_g/A_g;
% Yearly total electrical output in kWh/m2
total_electrical_output = etaE*total_I_sun*A_g/A_g;

pvt_collector_output.eta_tha=eta_tha;
pvt_collector_output.etaSum=etaSum;
pvt_collector_output.eta_E=eta_E;
pvt_collector_output.etaE=etaE;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.Tmean=Tmean;
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

f = figure(1023);
plot(xx,Qa,'DisplayName','Q_a')
ylabel('Thermal power output [W]');
yyaxis right
plot(xx,Qe,'DisplayName','Q_e')
ylabel('Electrical power output [W]');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Power Output Plot');
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
end
