function figures_cell_temperature_nonper(irradiance,cell_temp,ambient_temp,...
    WEATHER_output, TOOLBOX_input)
%FIGURES_CELL_TEMPERATURE_NONPER Generate cell temperature figures
% Show the results of the thermal module in case of non-periodic simulations.
%
% Parameters
% ----------
% irradiance : double
%   Irradiance for the period of interest
% cell_temp : cell
%   Temperature for each of the cells as a function of time
% ambient_temp : double
%   Ambient temperature for the period of interest
% WEATHER_output : struct
%   Simulation results of the WEATHER module
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by K Ganapathi Subramanian

indices = WEATHER_output.Period.dataset; 
mod_arr = TOOLBOX_input.Scene.Arrangement.ModuleArrangement;
[m,n] = size(mod_arr);
irr_to_plot = irradiance(indices(1):indices(2));
amb_temp = ambient_temp(indices(1):indices(2));
x = 1:length(irr_to_plot);

for i = 1:numel(mod_arr)
    temp_to_plot = cell_temp{i};
    temp_to_plot = temp_to_plot(indices(1):indices(2),:);
    temp_to_plot = mean(temp_to_plot,2);
    figure(36);
    subplot(m,n,i); 
    plot(x,temp_to_plot,'r-','LineWidth',2);
    hold on;
    plot(x,amb_temp','b-','LineWidth',2);
    hold off;
    yyaxis right
    fill([x, fliplr(x)],[irr_to_plot', zeros(1,length(irr_to_plot))],...
        'y','LineStyle','none');
    alpha(0.33);
    ylabel('Irradiance [W/m^2]');
    xlabel('Time [Hours]'); 
    ylabel('Temperature [^oC]');
    xlim([indices(1) indices(2)]);
    title("Module " + i);
    legend(gca,{'Mean Temperature',...
        'Ambient Temperature','Irradiance'});
end

for i = 1:numel(mod_arr)
    temp_to_plot = cell_temp{i};
    temp_to_plot = temp_to_plot(indices(1):indices(2),:);
    max_cell_temp = max(temp_to_plot,[],2);
    min_cell_temp = min(temp_to_plot,[],2);
    figure(37);
    subplot(m,n,i);
    fill([x, fliplr(x)],[max_cell_temp;flipud(min_cell_temp)],...
        [0.8 0.8 0.8],'LineStyle','none');
    hold on;
    plot(x,mean(temp_to_plot,2),'k-','LineWidth',2);
    hold off;
    xlabel('Time [Hours]'); 
    ylabel('Temperature [^oC]');
    xlim([indices(1) indices(2)]);
    title("Module " + i);
    legend(gca,{'Temperature Range','Mean Temperature'});
end
end
