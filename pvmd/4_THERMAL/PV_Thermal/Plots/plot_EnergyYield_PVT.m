%FIGURES_PVT_Output
% Parameters
% ----------
% irradiance : double
%   Irradiance for the period of interest
% day : double
%   Day of the month for each weather value
% month : double
%   Month for each weather value
% Developed by by ZUA.

function plot_EnergyYield_PVT(pvt_collector_output,WEATHER_output)

period=WEATHER_output.Period;
month=WEATHER_output.month;
day=WEATHER_output.day;
Irr=mean(WEATHER_output.Irr,2);
ambient_temperature=WEATHER_output.ambient_temperature;

Two = pvt_collector_output.Two;
Tc = pvt_collector_output.T;
% % % eta_tha = pvt_collector_output.eta_tha;
% % % eta_E = pvt_collector_output.eta_E;

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
    hours_day = 24;
    x = 1:time/hours_day;
    xx = 1:time;

% Xticks settings (per hour data)
    months = 24*[31,28,31,30,31,30,31,31,30,31,30,31];  
    ticks = 24*[15.5 31 45 59 74.5 90 105 120 135.5 151 166 181 197 212 227.5 ...
        243 258 273 288 303 317 332 350]-(period(1)+sum(months(1:period(2)-1)));
    month_labelss = {'Jan','','Feb','','Mar','','Apr','','May','','Jun',...
            '','Jul','','Aug','','Sep','','Oct','','Nov','','Dec'};
% Xticks settings (per day data)
    months = [31,28,31,30,31,30,31,31,30,31,30,31];  
    tick = [15.5 31 45 59 74.5 90 105 120 135.5 151 166 181 197 212 227.5 ...
        243 258 273 288 303 317 332 350]-(period(1)+sum(months(1:period(2)-1)));
    month_labels = {'Jan','','Feb','','Mar','','Apr','','May','','Jun',...
            '','Jul','','Aug','','Sep','','Oct','','Nov','','Dec'};

% Daily irradiance
    number_hours = length(Irr);
    number_days = number_hours/hours_day;
    daily_irradiance = mean(reshape(Irr,[hours_day,number_days]));

    [~,I] = max(daily_irradiance);
    highest_irradiance = Irr((I-1)*24+1:I*24);

% Highest irradiance day plot
    f = figure(32);
    area(1:24,highest_irradiance,'Linestyle','none','FaceColor','y');
    alpha(0.33);
    ylabel('Irradiance [W/m^2]'); xlabel('Time [Hours]');
    xlim([0 24])
    title(['Sunniest Day of the Period (',num2str(day(I*24)),'/',...
        num2str(month(I*24)),')']);
    SetFigureDefaults(f)

%% Daily irradiance plot
    f = figure(33);
    plot(x,daily_irradiance)
    if time/24 >=31
        xticks(tick); xticklabels(month_labels)
    else
        xlabel('Time [Days]')
    end
    ylabel('Daily Irradiation [kWh/m^2]')
    title('Daily Irradiation Plot');
    xlim([1 time/24])
    SetFigureDefaults(f)

end

%% Daily temperature plot
    f = figure(35);
    plot(xx,Two)
    hold on
    plot(xx,Tc)
    plot(xx,ambient_temperature)
    if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
    end
    ylabel('Temperature [^oC]')
    title('Daily Temperature Plot');
    xlim([1 time])
    legend('T_{fluid}', 'T_{cell}', 'T_{ambient}');
    SetFigureDefaults(f)

% % % % Daily efficiency plot
% % %     f = figure(36);
% % %     plot(xx,eta_tha)
% % %     ylabel('Thermal efficiency');
% % %     yyaxis right
% % %     plot(xx,eta_E)
% % %     if time/24 >=31
% % %         xticks(ticks); xticklabels(month_labelss)
% % %     else
% % %         xlabel('Time [hours]')
% % %     end
% % %     ylabel('Electrical efficiency');
% % %     title('Daily Efficiency Plot');
% % %     xlim([1 time])
% % %     legend('\eta_{thermal}', '\eta_{electrical}');
% % %     SetFigureDefaults(f)

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