function plot_pow_opt_graphs(Pdc_mod,Pdc_PO,eff_PO,eff_curve,Vmpp,...
    Pac,Pdc_system,model_PO,model_INV,periodic)
% Plot power optimizer graphs once the conversion module has finished
%
% Parameters
% ----------
% Pdc_mod : double/cell
%   DC power at module level
% Pdc_PO : double/cell
%   DC power at module level after the power optimizer
% eff_PO : double/cell
%   Efficiency of the power optimizer
% eff_curve : cell
%   Curve efficiencies characteristic of the power optimizer
% Vmpp : double
%   Voltage characteristics of the power optimizer
% Pac : double
%   AC power resulting from the central inverter
% Pdc_system : double/cell
%   DC power arriving to the central inverter
% model_PO : char
%   Power optimizer model
% model_INV : char
%   Central inverter model
% periodic : logical
%   Whether the simulations are periodic or not

% Power optimizer input and output
figure(); grid on
if periodic
    time = 1:size(Pdc_mod);
    bar(time,[Pdc_mod Pdc_PO],'FaceColor','flat','BarWidth', 1);
    title(['Power input and output of ' model_PO])
    xlabel('Time [h]'); ylabel('DC Power [W]');
    set(gca,'Fontsize',14);
    set(gcf,'Position',[200,200,800,400]);
else
    num_mods = length(Pdc_mod);
    time = 1:size(Pdc_mod{1},1);
    for i = 1:num_mods
        subplot(num_mods,1,i);
        bar(time,[Pdc_mod{i} Pdc_PO{i}],'FaceColor','flat');
        subplot_title = sprintf('Module %d',i);
        title(subplot_title);
        xlabel('Time [h]'); ylabel('DC Power [W]');
        set(gca,'Fontsize',12);
        set(gcf,'Position',[200,200,800,600]);
    end
end
legend('Power Optimizer Input Power','Power Optimizer Output Power')

% Time efficiency of the power optimizer
figure(); grid on;
title(['Efficiency of ' model_PO])
if periodic
    plot(time,eff_PO,'b')
    ylabel('Efficiency [%]')
    xlabel('Time [h]')
    set(gca,'Fontsize',14);
else
    for i = 1:num_mods
        subplot(num_mods,1,i);
        plot(time,eff_PO{i},'b')
        subplot_title = sprintf('Module %d',i);
        title(subplot_title);
        ylabel('Efficiency [%]'); xlabel('Time [h]')
        set(gca,'Fontsize',12);
    end
end
set(gcf, 'Position', [200 200 1200, 600]);

% Power optimize efficiency curves
figure(); hold on; grid minor; colors = 'brg';
legend_text = cell(3,1);
for i = 1:3
    plot(eff_curve{i}(:,1),eff_curve{i}(:,2),colors(i))
    legend_text{i} = sprintf('V_{mpp}^{%d} = %d V',i,Vmpp(i));
end
title(['Efficiency Curves of ' model_PO])
ylabel('Efficiency [%]'); xlabel('Input Power [W]')
legend(legend_text,'location','southeast')
set(gca,'Fontsize',14);

% Inverter input and output. Plot only is central inverter is selected
if ~isnan(model_INV)
    figure(); grid on;
    bar(1:length(Pdc_system),[Pdc_system Pac],'FaceColor','flat');
    legend('Inverter Input Power [V DC]','Inverter Output Power [V AC]')
    title(['Power input and output of ' model_INV])
    xlabel('Time [h]');
    ylabel('Po