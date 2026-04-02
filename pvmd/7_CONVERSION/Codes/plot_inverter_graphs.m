function plot_inverter_graphs(Pac, Pdc, eff, inv_type, periodic)
% Plot graphs once the conversion module has finished
%
% Parameters
% ----------
% Pac : double/cell
%   AC power resulting from the inverter
% Pdc : double/cell
%   DC power
% eff : double/cell
%   Efficiency of the inverter
% periodic : logical
%   Whether the simulations are periodic or not

% Inverter input and output power
figure(); grid on;
if periodic || strcmp(inv_type,'CEN')
    P_length = size(Pdc,1);
    time = 1:P_length;
    bar(time,[Pdc Pac],'FaceColor','flat','BarWidth', 1);
    ylim([0 max(Pdc+5)]); xlim([0 P_length]);
    xlabel('Time [h]'); ylabel('Output [W]');
    set(gca,'Fontsize',14);
    set(gcf,'Position',[200,200,800,400]);
else %non-periodic simulations for micro/string inverters
    time = 1:size(Pdc{1},1);
    parallel = length(Pdc);
    if strcmp(inv_type,'MIC')
        block = 'module';
    elseif strcmp(inv_type,'STR')
        block = 'string';
    end
    for i = 1:parallel
        subplot(parallel,1,i);
        bar(time,[Pdc{i} Pac{i}],'FaceColor','flat');
        ylim([0 max(Pdc{i}+5)]);
        subplot_title = sprintf('Inverter performance %s %d',block,i);
        title(subplot_title);
        xlabel('Time [h]'); ylabel('Output [W]');
        set(gca,'Fontsize',12);
    end
    set(gcf,'Position',[200,200,800,600]);
end
legend(gca,{'Inverter Input Power','Inverter Output Power'});


% Inverter efficiency
figure(); grid on;
if periodic || strcmp(inv_type,'CEN')
    plot(time,eff,'b');
    xlim([0 P_length]);
    ylabel('Efficiency [%]'); xlabel('Time [h]');
    set(gca,'Fontsize',14);
else
    for i = 1:parallel
        subplot(parallel,1,i);
        plot(time,eff{i},'b');
        subplot_title = sprintf('Efficiency of inverter %s %d',block,i);
        title(subplot_title);
        ylabel('Efficiency [%]'); xlabel('Time [h]');
        set(gca,'Fontsize',12);