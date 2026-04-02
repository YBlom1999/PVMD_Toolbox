function figures_cell_temperature_per(irradiance,cell_temp,ambient_temp,day,month,period)
%FIGURES_CELL_TEMPERATURE_PER Generate cell temperature figures
% Show the results of the thermal module in case of periodic simulations.
%
% Parameters
% ----------
% irradiance : double
%   Irradiance for the period of interest
% cell_temp : double
%   Temperature for each of the cells as a function of time
% ambient_temp : double
%   Ambient temperature for the period of interest
% day : double
%   Day of the month for each weather value
% month : double
%   Month for each weather value
%
% Developed by unknown (A. Nour?). Commented by A. Alcaniz

irradiance_avg = mean(irradiance,2); %Calculate the average irradiance over the module
time = length(irradiance_avg);
if time <= 168
    x = 1:time;
    
    % Cell temperature plot
    f = figure(31); hold on;
    max_cell_temp = max(cell_temp,[],2); 
    min_cell_temp = min(cell_temp,[],2);
    fill([x, fliplr(x)],[max_cell_temp;flipud(min_cell_temp)],...
        [0.8 0.8 0.8],'LineStyle','none');
    plot(x,mean(cell_temp,2),'k','LineWidth',2);
    hold off
    
    xlabel('Time [Hours]'); ylabel('Temperature [^oC]')
    legend(gca,{'Cell Temperature','Average Temperature'})
    title('Cell temperatures');
    xlim([1 time])
    SetFigureDefaults(f)
    
    % Irradiance and temperature plots
    f = figure(33); hold on;
    
    yyaxis right
    fill([x, fliplr(x)],[irradiance_avg; zeros(time,1)],'y','LineStyle','none')
    alpha(0.33)
    ylabel('Irradiance [W/m^2]');
    f.CurrentAxes.YColor = 'k';
    
    yyaxis left
    plot(x,mean(cell_temp,2),'r-','LineWidth',2);
    plot(x,ambient_temp,'b-','LineWidth',2);
    f.CurrentAxes.YColor = 'k';
    hold off
    
    xlabel('Time [Hours]'); ylabel('Temperature [^oC]');
    legend({'T_{module}','T_{amb}','Irradiance'},'Location','NorthEast',...
        'FontSize',14);
    xlim([1 time])
    SetFigureDefaults(f)
else
    hours_day = 24;
    x = 1:time/hours_day;
    
    % Mean, min and max daily temperatures
    [number_hours,number_cells] = size(cell_temp);
    number_days = number_hours/hours_day; N = number_cells*hours_day;
    daily_cell_temp = reshape(cell_temp',[N,number_days]);
    mean_daily_cell_temp = mean(daily_cell_temp);
    max_daily_cell_temp = max(daily_cell_temp);
    min_daily_cell_temp = min(daily_cell_temp);
    
    % Xticks settings
    months = [31,28,31,30,31,30,31,31,30,31,30,31];  
    tick = [15.5 31 45 59 74.5 90 105 120 135.5 151 166 181 197 212 227.5 ...
        243 258 273 288 303 317 332 350]-(period(1)+sum(months(1:period(2)-1)));
    month_labels = {'Jan','','Feb','','Mar','','Apr','','May','','Jun',...
            '','Jul','','Aug','','Sep','','Oct','','Nov','','Dec'};
    
    % Daily irradiance
    daily_irradiance = mean(reshape(irradiance_avg,[hours_day,number_days]));
    [~,I] = max(daily_irradiance);
    highest_irradiance = irradiance_avg((I-1)*24+1:I*24);
    
    
    % Module temperature plot
    f = figure(31); hold on;
    fill([x, fliplr(x)],[min_daily_cell_temp, ...
        fliplr(max_daily_cell_temp)],[0.8 0.8 0.8],'LineStyle','none')
    plot(x,mean_daily_cell_temp,'k','LineWidth',1.5)
    hold off
    
    if time/24>=31
        xticks(tick); xticklabels(month_labels)
    else
        xlabel('Time [Days]')
    end
    ylabel('Temperature [^oC]'); xlim([1 time/24])
    title('Module Temperature');
    legend(gca,{'Cell Temperature','Average Temperature'})
    SetFigureDefaults(f)
    
    
    % Highest irradiance day plot
    f = figure(32);
    area(1:24,highest_irradiance,'Linestyle','none','FaceColor','y');
    alpha(0.33);
    
    ylabel('Irradiance [W/m^2]'); xlabel('Time [Hours]');
    xlim([0 24])
    title(['Sunniest Day of the Period (',num2str(day(I*24)),'/',...
        num2str(month(I*24)),')']);
    SetFigureDefaults(f)
    
    
    % Daily irradiance plot
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

function SetFigureDefaults(f)
    % Set several figure characteristics for correct representation
    %
    % INPUT
    % f    figure object to be modified
    
    f.CurrentAxes.FontSize = 14;
    f.CurrentAxes.Box = 'on';
    f.CurrentAxes.XGrid = 'on';
    f.CurrentAxes.YGrid = 'on';
    f.CurrentAxes.XLabel.FontSize = 16;
    f.CurrentAxes.YLabel.FontSize = 16;
    f.CurrentAxes.Title.FontSize = 16;
end

end

