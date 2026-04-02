function [RMSE_J, NRMSE_J, RMSE_V,NRMSE_V ]= plot_IVcurveDifference(V1,I1,V2,I2,Curve_names,name)

%V in volt
%I in A/m2
%expects V1,I1 to be simulated IV curve
%expects V2,I2 to be measured IV curve

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


plot(V1,I1);
plot(V2,I2,'o');

%Legends - you can choose which plots appear in the legend
lgd=legend(my_ax1,Curve_names,'Location','southwest');

%Adjusting axes
box(my_ax1,'on'); %draw a box around the axes
my_ax1.LineWidth = 1;%change the thickness of the box and axes ticks
%axis(my_ax1,'equal'); %use the same spacing for x and y grids

my_ax1.XLim = [0 max(max(V1),max(V2))*1.05]; %change axes limits
my_ax1.YLim = [0 max(max(I1),max(I2))*1.05]; %change axes limits
%Change the ticks displayed on the x and y axes
%my_ax1.XTick = [-1 0 1];
%my_ax1.YTick = [-1 -0.5 0 0.5 1];
%my_ax1.XTickLabel = {'-1.0','0.0','+1.0'};%you can replace the tick labels

my_ax1.XLabel.String = 'Voltage [V]';
my_ax1.YLabel.String = 'Current [A/m2]';
my_ax1.XLabel.FontWeight = 'Bold';
my_ax1.YLabel.FontWeight = 'Bold';

my_ax1.FontName = 'Calibri';
my_ax1.FontSize = 14;

%% prepare for deviation calculation

V_end=min(max(V1),max(V2));
V_end_2digts=round(V_end,2)-0.01;%make sure it is rounded down(prevents Nan in interpolation)
 
V_intpl=0:0.01:V_end_2digts;

ind_Vpositive=find(V1>0);

Jsc_sim_intpl=interp1(V1(ind_Vpositive),I1(ind_Vpositive),V_intpl);
Jsc_meas_intpl=interp1(V2,I2,V_intpl);

ind_not_nan=find(~isnan(Jsc_sim_intpl));
Jsc_sim_intpl=Jsc_sim_intpl(ind_not_nan);
Jsc_meas_intpl=Jsc_meas_intpl(ind_not_nan);

%AbsMeanDev_J = mean(abs(Jsc_sim_intpl - Jsc_meas_intpl));
%RelMeanDev_J = mean(abs(Jsc_sim_intpl - Jsc_meas_intpl)./Jsc_meas_intpl);
RMSE_J= sqrt(sum((Jsc_sim_intpl-Jsc_meas_intpl).^2)/length(Jsc_sim_intpl));
NRMSE_J= sqrt(sum((Jsc_sim_intpl-Jsc_meas_intpl).^2)/(length(Jsc_sim_intpl)*mean(Jsc_meas_intpl)));

J_end=min(max(I1),max(I2));
J_end_2digts=round(J_end,2)-0.01;%make sure it is rounded down(prevents Nan in interpolation)
J_intpl=0:0.01:J_end_2digts;

V_sim_intpl=interp1(I1,V1,J_intpl);
[~, ind_unique, ~]=unique(I2);
V_meas_intpl=interp1(I2(ind_unique),V2(ind_unique),J_intpl);

%AbsMeanDev_V = mean(abs(V_sim_intpl - V_meas_intpl));
%RelMeanDev_V = mean(abs(V_sim_intpl - V_meas_intpl)./V_meas_intpl);
RMSE_V= sqrt(sum((V_sim_intpl-V_meas_intpl).^2)/length(V_sim_intpl));
NRMSE_V= sqrt(sum((V_sim_intpl-V_meas_intpl).^2)/(length(V_sim_intpl)*mean(V_meas_intpl)));

%Legends - you can choose which plots appear in the legend

%title(lgd,'Relative deviation'+string(RelMeanDev_J));

savefig(name+'.fig');
close all;