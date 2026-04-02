function [my_fig] = plot_losses(Plot_Components,Plot_Power,Plot_Categories,Plot_Power_Categories,Plot_type,title)
%plot_losses Creates a figure with the results
%
% This function plots the results of the simulation.
%
% Parameters
% ----------
% Plot_Components : Char
%   The names of all the components
% Plot_Power : Char
%   The values of all the components
% Plot_Categories : Char
%   The names of all the categories
% Plot_Power_Categories : Char
%   The values of all the categories
% Plot_type : Char
%   The type of simulation (DC or AC/ single junction or tandem)
% title : Char
%   The title of the figure
%
% Returns
% -------
% my_fig : figure
%   The figure of the results.
%
% Developed by Y. Blom

my_fig = figure();
sgtitle(title,'Fontsize',20)
fig_width_cm = 15;%width of the figure in cm
fig_height_cm = 18;%width of the figure in cm
fig_width = round(fig_width_cm/(2.54/96)); %in pixels
fig_height= round(fig_height_cm/(2.54/96)); %in pixels
x_start = 200;
y_start = 50;
my_fig.Position = [x_start y_start fig_width fig_height];
% Line colors in normalized RGB coordinates
color1 = [0.8500    0.3250    0.0980];
color2 = [0    0.4470    0.7410];
color3 = [0.9290    0.6940    0.1250];
color4 = [0.4940    0.1840    0.5560];
color5 = [0.4660    0.6740    0.1880];

result_plot = subplot(1,1,1);
result_plot.Position = result_plot.Position + [0.25 -0.04 -0.25 0];
result_plot.Clipping = 'off';
hold(result_plot,'on')
box(result_plot,'on')

yb = (1:length(Plot_Power));
xb = Plot_Power;

%% Situation AC single
if Plot_type == 1
    bplot_01 = barh(result_plot,yb(1),xb(1),'FaceColor',color5,'EdgeAlpha',0);
    bplot_02 = barh(result_plot,yb(2),xb(2),'FaceColor',color4,'EdgeAlpha',0);
    bplot_03 = barh(result_plot,yb(3),xb(3),'FaceColor',color4,'EdgeAlpha',0);
    bplot_04 = barh(result_plot,yb(4),xb(4),'FaceColor',color4,'EdgeAlpha',0);
    bplot_05 = barh(result_plot,yb(5),xb(5),'FaceColor',color4,'EdgeAlpha',0);
    bplot_06 = barh(result_plot,yb(6),xb(6),'FaceColor',color3,'EdgeAlpha',0);
    bplot_07 = barh(result_plot,yb(7),xb(7),'FaceColor',color3,'EdgeAlpha',0);
    bplot_08 = barh(result_plot,yb(8),xb(8),'FaceColor',color3,'EdgeAlpha',0);
    bplot_09 = barh(result_plot,yb(9),xb(9),'FaceColor',color3,'EdgeAlpha',0);
    bplot_10 = barh(result_plot,yb(10),xb(10),'FaceColor',color2,'EdgeAlpha',0);
    bplot_11 = barh(result_plot,yb(11),xb(11),'FaceColor',color2,'EdgeAlpha',0);
    bplot_12 = barh(result_plot,yb(12),xb(12),'FaceColor',color2,'EdgeAlpha',0);
    bplot_13 = barh(result_plot,yb(13),xb(13),'FaceColor',color2,'EdgeAlpha',0);
    bplot_14 = barh(result_plot,yb(14),xb(14),'FaceColor',color1,'EdgeAlpha',0);
    bplot_15 = barh(result_plot,yb(15),xb(15),'FaceColor',color1,'EdgeAlpha',0);
    bplot_16 = barh(result_plot,yb(16),xb(16),'FaceColor',color1,'EdgeAlpha',0);
    bplot_17 = barh(result_plot,yb(17),xb(17),'FaceColor',color1,'EdgeAlpha',0);
    bplot_18 = barh(result_plot,yb(18),xb(18),'FaceColor',color1,'EdgeAlpha',0);
    bplot_19 = barh(result_plot,yb(19),xb(19),'FaceColor',color1,'EdgeAlpha',0);


    ylim(result_plot,[0,length(Plot_Power)+1])
    xlabel(result_plot, 'Percentage [%]','FontSize',11)
    result_plot.YTick = [1:length(Plot_Power)];
    result_plot.YTickLabel = categorical(Plot_Components);
    result_plot.YTickLabelRotation = 0;
    result_plot.FontSize = 15;
    text_result_per = [num2str(xb,'%0.1f'),char(ones(length(xb),1)*37)];
    text_result_groups = [num2str(Plot_Power_Categories,'%0.1f'),char(ones(length(Plot_Power_Categories),1)*37)];

    result_legend_text = char(Plot_Categories');
    result_legend = legend(result_plot,[bplot_19 bplot_13 bplot_09 bplot_05 bplot_01],{result_legend_text}, 'Location','southoutside', 'box','off','FontSize',15);
    text(result_plot,max(xb,0)+5,yb,text_result_per,'vert','middle','horiz','center','FontSize',15,'FontName', 'Calibri','Rotation' ,0);
    if min(xb) < -10
        xlim(result_plot,[-12, 45])
        rectangle(result_plot,'position',[-7 -12 57 8],'LineWidth',1.5)
        text(result_plot,44,-8,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    elseif min(xb) < -7
        xlim(result_plot,[-10, 45])
        rectangle(result_plot,'position',[-4 -12 54 8],'LineWidth',1.5)
        text(result_plot,44,-8,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    elseif min(xb) < -5
        xlim(result_plot,[-7, 45])
        rectangle(result_plot,'position',[-2 -12 51 8],'LineWidth',1.5)
        text(result_plot,44,-8,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    elseif min(xb) < -2
        xlim(result_plot,[-5, 45])
        rectangle(result_plot,'position',[-2 -12 51 8],'LineWidth',1.5)
        text(result_plot,44,-8,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    else
        xlim(result_plot,[-2, 50])
        rectangle(result_plot,'position',[-1 -12 56 8],'LineWidth',1.5)
        text(result_plot,50,-8,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')

    end
end
%% Situation AC tandem
if Plot_type == 2
    bplot_01 = barh(result_plot,yb(1),xb(1),'FaceColor',color5,'EdgeAlpha',0);
    bplot_02 = barh(result_plot,yb(2),xb(2),'FaceColor',color5,'EdgeAlpha',0);
    bplot_03 = barh(result_plot,yb(3),xb(3),'FaceColor',color4,'EdgeAlpha',0);
    bplot_04 = barh(result_plot,yb(4),xb(4),'FaceColor',color4,'EdgeAlpha',0);
    bplot_05 = barh(result_plot,yb(5),xb(5),'FaceColor',color4,'EdgeAlpha',0);
    bplot_06 = barh(result_plot,yb(6),xb(6),'FaceColor',color4,'EdgeAlpha',0);
    bplot_07 = barh(result_plot,yb(7),xb(7),'FaceColor',color3,'EdgeAlpha',0);
    bplot_08 = barh(result_plot,yb(8),xb(8),'FaceColor',color3,'EdgeAlpha',0);
    bplot_09 = barh(result_plot,yb(9),xb(9),'FaceColor',color3,'EdgeAlpha',0);
    bplot_10 = barh(result_plot,yb(10),xb(10),'FaceColor',color3,'EdgeAlpha',0);
    bplot_11 = barh(result_plot,yb(11),xb(11),'FaceColor',color2,'EdgeAlpha',0);
    bplot_12 = barh(result_plot,yb(12),xb(12),'FaceColor',color2,'EdgeAlpha',0);
    bplot_13 = barh(result_plot,yb(13),xb(13),'FaceColor',color2,'EdgeAlpha',0);
    bplot_14 = barh(result_plot,yb(14),xb(14),'FaceColor',color2,'EdgeAlpha',0);
    bplot_15 = barh(result_plot,yb(15),xb(15),'FaceColor',color1,'EdgeAlpha',0);
    bplot_16 = barh(result_plot,yb(16),xb(16),'FaceColor',color1,'EdgeAlpha',0);
    bplot_17 = barh(result_plot,yb(17),xb(17),'FaceColor',color1,'EdgeAlpha',0);
    bplot_18 = barh(result_plot,yb(18),xb(18),'FaceColor',color1,'EdgeAlpha',0);
    bplot_19 = barh(result_plot,yb(19),xb(19),'FaceColor',color1,'EdgeAlpha',0);
    bplot_20 = barh(result_plot,yb(20),xb(20),'FaceColor',color1,'EdgeAlpha',0);


    ylim(result_plot,[0,length(Plot_Power)+1])
    xlabel(result_plot, 'Percentage [%]','FontSize',11)
    result_plot.YTick = [1:length(Plot_Power)];
    result_plot.YTickLabel = categorical(Plot_Components);
    result_plot.YTickLabelRotation = 0;
    result_plot.FontSize = 15;
    text_result_per = [num2str(xb,'%0.1f'),char(ones(length(xb),1)*37)];
    text_result_groups = [num2str(Plot_Power_Categories,'%0.1f'),char(ones(length(Plot_Power_Categories),1)*37)];

    result_legend_text = char(Plot_Categories');
    result_legend = legend(result_plot,[bplot_20 bplot_14 bplot_10 bplot_06 bplot_01],{result_legend_text}, 'Location','southoutside', 'box','off','FontSize',15);
    text(result_plot,max(xb,0)+5,yb,text_result_per,'vert','middle','horiz','center','FontSize',15,'FontName', 'Calibri','Rotation' ,0);
    if min(xb) < -12
        xlim(result_plot,[-15, 40])
        rectangle(result_plot,'position',[-9 -12.5 55 8],'LineWidth',1.5)
        text(result_plot,40,-8.5,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    elseif min(xb) < -10
        xlim(result_plot,[-12, 40])
        rectangle(result_plot,'position',[-7 -12.5 54 8],'LineWidth',1.5)
        text(result_plot,40,-8.5,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    elseif min(xb) < -7
        xlim(result_plot,[-10, 40])
        rectangle(result_plot,'position',[-5 -12.5 50 8],'LineWidth',1.5)
        text(result_plot,40,-8.5,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    elseif min(xb) < -5
        xlim(result_plot,[-7, 40])
        rectangle(result_plot,'position',[-3 -12.5 48 8],'LineWidth',1.5)
        text(result_plot,40,-8.5,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    elseif min(xb) < -2
        xlim(result_plot,[-5, 40])
        rectangle(result_plot,'position',[-1 -12.5 46 8],'LineWidth',1.5)
        text(result_plot,40,-8.5,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    else
        xlim(result_plot,[-2, 50])
        rectangle(result_plot,'position',[3 -12.5 53 8],'LineWidth',1.5)
        text(result_plot,50,-8.5,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')

    end
end
%% Situation DC single
if Plot_type == 3
    bplot_01 = barh(result_plot,yb(1),xb(1),'FaceColor',color5,'EdgeAlpha',0);
    bplot_02 = barh(result_plot,yb(2),xb(2),'FaceColor',color4,'EdgeAlpha',0);
    bplot_03 = barh(result_plot,yb(3),xb(3),'FaceColor',color4,'EdgeAlpha',0);
    bplot_04 = barh(result_plot,yb(4),xb(4),'FaceColor',color3,'EdgeAlpha',0);
    bplot_05 = barh(result_plot,yb(5),xb(5),'FaceColor',color3,'EdgeAlpha',0);
    bplot_06 = barh(result_plot,yb(6),xb(6),'FaceColor',color3,'EdgeAlpha',0);
    bplot_07 = barh(result_plot,yb(7),xb(7),'FaceColor',color3,'EdgeAlpha',0);
    bplot_08 = barh(result_plot,yb(8),xb(8),'FaceColor',color2,'EdgeAlpha',0);
    bplot_09 = barh(result_plot,yb(9),xb(9),'FaceColor',color2,'EdgeAlpha',0);
    bplot_10 = barh(result_plot,yb(10),xb(10),'FaceColor',color2,'EdgeAlpha',0);
    bplot_11 = barh(result_plot,yb(11),xb(11),'FaceColor',color2,'EdgeAlpha',0);
    bplot_12 = barh(result_plot,yb(12),xb(12),'FaceColor',color1,'EdgeAlpha',0);
    bplot_13 = barh(result_plot,yb(13),xb(13),'FaceColor',color1,'EdgeAlpha',0);
    bplot_14 = barh(result_plot,yb(14),xb(14),'FaceColor',color1,'EdgeAlpha',0);
    bplot_15 = barh(result_plot,yb(15),xb(15),'FaceColor',color1,'EdgeAlpha',0);
    bplot_16 = barh(result_plot,yb(16),xb(16),'FaceColor',color1,'EdgeAlpha',0);
    bplot_17 = barh(result_plot,yb(17),xb(17),'FaceColor',color1,'EdgeAlpha',0);

    ylim(result_plot,[0,length(Plot_Power)+1])
    xlabel(result_plot, 'Percentage [%]','FontSize',11)
    result_plot.YTick = [1:length(Plot_Power)];
    result_plot.YTickLabel = categorical(Plot_Components);
    result_plot.YTickLabelRotation = 0;
    result_plot.FontSize = 15;
    text_result_per = [num2str(xb,'%0.1f'),char(ones(length(xb),1)*37)];
    text_result_groups = [num2str(Plot_Power_Categories,'%0.1f'),char(ones(length(Plot_Power_Categories),1)*37)];

    result_legend_text = char(Plot_Categories');
    result_legend = legend(result_plot,[bplot_17 bplot_11 bplot_07 bplot_03 bplot_01],{result_legend_text}, 'Location','southoutside', 'box','off','FontSize',15);
    text(result_plot,max(xb,0)+5,yb,text_result_per,'vert','middle','horiz','center','FontSize',15,'FontName', 'Calibri','Rotation' ,0);
    if min(xb) < -2
        xlim(result_plot,[-5, 40])
        rectangle(result_plot,'position',[-1 -10.7 46 7],'LineWidth',1.5)
    else
        xlim(result_plot,[-2, 50])
        rectangle(result_plot,'position',[2.3 -10.7 53 7],'LineWidth',1.5)
        text(result_plot,50,-7,text_result_groups,'vert','middle','horiz','center','FontSize',15,'FontName', 'Helvetica')
    end

end
%% Situation DC tandem
if Plot_type == 4
    bplot_01 = barh(result_plot,yb(1),xb(1),'FaceColor',color5,'EdgeAlpha',0);
    bplot_02 = barh(result_plot,yb(2),xb(2),'FaceColor',color5,'EdgeAlpha',0);
    bplot_03 = barh(result_plot,yb(3),xb(3),'FaceColor',color4,'EdgeAlpha',0);
    bplot_04 = barh(result_plot,yb(4),xb(4),'FaceColor',color4,'EdgeAlpha',0);
    bplot_05 = barh(result_plot,yb(5),xb(5),'FaceColor',color3,'EdgeAlpha',0);
    bplot_06 = barh(result_plot,yb(6),xb(6),'FaceColor',color3,'EdgeAlpha',0);
    bplot_07 = barh(result_plot,yb(7),xb(7),'FaceColor',color3,'EdgeAlpha',0);
    bplot_08 = barh(result_plot,yb(8),xb(8),'FaceColor',color3,'EdgeAlpha',0);
    bplot_09 = barh(result_plot,yb(9),xb(9),'FaceColor',color2,'EdgeAlpha',0);
    bplot_10 = barh(result_plot,yb(10),xb(10),'FaceColor',color2,'EdgeAlpha',0);
    bplot_11 = barh(result_plot,yb(11),xb(11),'FaceColor',color2,'EdgeAlpha',0);
    bplot_12 = barh(result_plot,yb(12),xb(12),'FaceColor',color2,'EdgeAlpha',0);
    bplot_13 = barh(result_plot,yb(13),xb(13),'FaceColor',color1,'EdgeAlpha',0);
    bplot_14 = barh(result_plot,yb(14),xb(14),'FaceColor',color1,'EdgeAlpha',0);
    bplot_15 = barh(result_plot,yb(15),xb(15),'FaceColor',color1,'EdgeAlpha',0);
    bplot_16 = barh(result_plot,yb(16),xb(16),'FaceColor',color1,'EdgeAlpha',0);
    bplot_17 = barh(result_plot,yb(17),xb(17),'FaceColor',color1,'EdgeAlpha',0);
    bplot_18 = barh(result_plot,yb(18),xb(18),'FaceColor',color1,'EdgeAlpha',0);

    ylim(result_plot,[0,length(Plot_Power)+1])
    xlabel(result_plot, 'Percentage [%]','FontSize',11)
    result_plot.YTick = [1:length(Plot_Power)];
    result_plot.YTickLabel = categorical(Plot_Components);
    result_plot.YTickLabelRotation = 0;
    result_plot.FontSize = 15;
    text_result_per = [num2str(xb,'%0.1f'),char(ones(length(xb),1)*37)];
    text_result_groups = [num2str(Plot_Power_Categories,'%0.1f'),char(ones(length(Plot_Power_Categories),1)*37)];

    result_legend_text = char(Plot_Categories');
    result_legend = legend(result_plot,[bplot_18 bplot_12 bplot_08 bplot_04 bplot_01],{result_legend_text}, 'Location','southoutside', 'box','off','FontSize',15);
    text(result_plot,max(xb,0)+5,yb,text_result_per,'vert','middle','horiz','center','FontSize',15,'FontName', 'Calibri','Rotation' ,0);
    if min(xb)< -2
        xlim(result_plot,[-5, 40])
        rectangle(result_plot,'position',[11 -11.3 18 7.5],'LineWidth',1.5)
    e