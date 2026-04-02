function [pvt_collector_output] = Solar_PV_thermal_II(~, WEATHER_output)
% Calculates the thermal and PV yield of PVT collector with unique coefficients
% Thermo photovoltaic collector - pap-II - 1.63 m2
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

A_g = 1.6269;                   % Gross Collector Area (m2)
dm_f = 0.027;                   % Fluid Mass Flowrate (kg/s)

% Second Order Thermal Efficiency Equation unique coefficients
a1a = 4.2755;                   % unique coefficient W/(m2K)                                          
etao_tha = 0.5447;              % unique coefficient
C_f = 4186;                     % specific heat (J/kg°C) 
% Electrical Efficiency Equation unique coefficients
e1a = 1.1618;                   % unique coefficient W/(m2K) 
etao_ele = 0.1525;              % unique coefficient

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
NOCT = 45;                       % typical module at 48°C (best module operated at a NOCT of 33°C, the worst at 58°C)
S = I*0.1;                       % insolation in mW/cm^2
Tc = Tam + S*(NOCT - 20)/80;
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
