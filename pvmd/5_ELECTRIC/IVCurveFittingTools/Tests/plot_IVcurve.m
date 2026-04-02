function plot_IVcurve(V,I,Curve_names,name)

%create a new figure and a new pair of axes
my_fig1 = figure();
my_ax1 = axes(my_fig1);

%define the width and the H/W ratio according to your preferences
fig_width_cm = 12;%width of the figure in cm
fig_h_w_ratio = 1;

%modify the dimensions of the figure
fig_width = round(fig_width_cm/(2.54/96)); %in pixels
fig_height= round(fig_width/fig_h_w_ratio);
x_start = 400;
y_start = 400;
my_fig1.Position = [x_start y_start fig_width fig_height];

hold(my_ax1,'on');%to retain current plot when adding new plots


for i=1:length(V(:,1))
    plot(V(i,:),I);
end

%Legends - you can choose which plots appear in the legend
lgd=legend(my_ax1,Curve_names,'Location','southwest');

%Adjusting axes
box(my_ax1,'on'); %draw a box around the axes
my_ax1.LineWidth = 1;%change the thickness of the box and axes ticks
%axis(my_ax1,'equal'); %use the same spacing for x and y grids
my_ax1.XLim = [0 max(max(V))*1.05]; %change axes limits
my_ax1.YLim = [0 max(I)*1.05]; %change axes limits
%Change the ticks displayed on the x and y axes
%my_ax1.XTick = [-1 0 1];
%my_ax1.YTick = [-1 -0.5 0 0.5 1];
%my_ax1.XTickLabel = {'-1.0','0.0','+1.0'};%you can replace the tick labels

my_ax1.XLabel.String = 'Voltage [V]';
my_ax1.YLabel.String = 'Current [A]';
my_ax1.XLabel.FontWeight = 'Bold';
my_ax1.YLabel.FontWeight = 'Bold';

my_ax1.FontName = 'Calibri';
my_ax1.FontSize = 14;


savefig(name+'.fig');
close all;