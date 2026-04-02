function plot_IVCurve(period,V_module,I,TOOLBOX_input)
% Plots the IV of the selected hour
%
%---INPUT---
% period - simulate time period typically WEATHER_output.Period
% V_module - Module voltage as functon of current
% I - Interval of current to be plotted
% TOOLBOX_input - the input of the user

%---OUTPUT---
%(none)
Months=[31,28,31,30,31,30,31,31,30,31,30,31];
Q = {'Month','Day','Hour'};
default = {'1','1','12'};
options.Resize='on';
A = str2double(inputdlg(Q,'Detailed Plots',1,default,options));
cell1=sum(Months(1:period(2)-1))*24+(period(1)-1)*24;
cell2=sum(Months(1:period(4)-1))*24+(period(3))*24;

req=sum(Months(1:A(1)-1))*24+(A(2)-1)*24+A(3);
if req < cell2 && req > cell1
%     GEO=str2double(evalin('base','MODULE_output.geometry'));
    cellnum=req-cell1;
%     if TYPE == 'T-F'
%     temperature=T_(cellnum,1:end)';
%     CellCurrent=J_(cellnum,1:end)';
%     else
%     temperature=reshape(T_(cellnum,1:end),GEO(1),GEO(2))';
%     for i=1:length(J_(1,1,:))
%     CellCurrent=Acell*reshape(J_(cellnum,1:end,i),GEO(1),GEO(2));
%     h=heatmap(CellCurrent(1,:),'Title','Cell Photo-generated Current [A]')
% Ax=gca
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% colormap hot
% set(gcf,'Position',[300 300 800 380])
% mycmap = get(gca,'Colormap');
% caxis([0 10])
%     end  
%     
% figure(220)
% clo(figure(220))
% h=heatmap(temperature-273.15,'Title','Cell Temperature [?C]')
% Ax=gca
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% colormap jet
% set(gcf,'Position',[300 300 800 380])
% mycmap = get(gca,'Colormap');
% caxis([-10 70])
% figure(221)
% clo(figure(221))
% if Model==1
%     CellCurrent(1,:)=CellCurrent(1,:)*Isc_STC/(408*Acell);
% end
% h=heatmap(CellCurrent(1,:),'Title','Cell Photo-generated Current [A]')
% Ax=gca
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% colormap hot
% set(gcf,'Position',[300 300 800 380])
% mycmap = get(gca,'Colormap');
% caxis([0 10])
time=cell2-cell1;
figure(222)
if TOOLBOX_input.runPeriodic
    if strcmp(TOOLBOX_input.electric.TYPE,'Non-REC')
        plot(V_module(cellnum,:),I);
    else
        plot(V_module(cellnum,:),I(cellnum,:));
    end
    hold on
    if length(V_module(:,1))>time
        if strcmp(TOOLBOX_input.electric.TYPE,'Non-REC')
            plot(V_module(cellnum,:),I);
        else
            plot(V_module(cellnum+time,:),I(cellnum+time,:));
        end
        legend('Top Cell','Bottom Cell')
    end
    box on
    xlabel('Voltage [V]')
    ylabel('Current [A]')
    hold off
    title("I-V Curve of the Module");
    grid on;
else %non-periodic simulations
    for i = 1:numel(I)
        v = V_module{i};
        i1 = I{i};
        subplot(size(I,1),size(I,2),i);
        plot(v(cellnum,:),i1);
            
        hold on;
        if length(v(:,1))>time
            if strcmp(MODULE_output.TYPE,'Non-REC')
                plot(v(cellnum,:),i1);
            else
                plot(v(cellnum+time,:),i1(cellnum+time,:));
            end
            legend(gca,{'Top Cell','Bottom Cell'})
        end
        box on
        xlabel('Voltage [V]')
        ylabel('Current [A]')
        hold off;
        title("I-V Curve for Module: " + i);
        grid on;
    end
end
end
