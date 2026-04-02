function [pvt_collector_output] = Solar_thermal_air(~, WEATHER_output)
% Calculates the thermal yield of thermal collector with unique coefficient
% Glazed Flat Plate Solar Air Heating Collector - 2.6 m2
% OG-100 Solar Thermal Collector Certification (10002148 - Trigo Energies Inc - SLG - Glazed Flat Plate)
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

A_g = 2.6;                      % Gross Collector Area (m2)
dm_f = 0.0844;                  % Fluid Mass Flowrate (kg/m2s)

% Second Order Thermal Efficiency Equation unique coefficients
a1a = 6.48;                     % W/(m2K)                                          
a2a = 0.028;                    % W/(m2K2) 
etao_tha = 0.40;
C_f = 1007;                     % specific heat (J/kg*C) 

I = mean(WEATHER_output.Irr,2)';
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
total_I_sun = sum(mean(WEATHER_output.Irr,2))/1000; 
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

f = figure(1022);
plot(xx,Qa,'DisplayName','Q_a')
ylabel('power output [W]');
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
