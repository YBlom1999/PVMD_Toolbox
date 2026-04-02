%create a new figure and a new pair of axes
my_fig1 = figure();
my_ax1 = axes(my_fig1);

%adapt this ratio according to your preferences
fig_width_cm = 18;%width of the figure in cm
fig_height_cm = 9;%width of the figure in cm
fig_width = round(fig_width_cm/(2.54/96)); %in pixels
fig_height= round(fig_height_cm/(2.54/96)); %in pixels
x_start = 800;
y_start = 400;
my_fig1.Position = [x_start y_start fig_width fig_height];

hold(my_ax1,'on');%to retain current plot when adding new plots

%numberOfBars=2;
Months=[31,28,31,30,31,30,31,31,30,31,30,31];

%monthlyValue=zeros(length(Months),numberOfBars);
monthsInHours=zeros(12,1);
monthsInHours(1)= 24*31/2;


for i=2:length(Months)
    monthsInHours(i)=monthsInHours(i-1)+24*(Months(i)+Months(i-1))/2;
end


mismatch=100*(Electric4T.Impp(:,1)-Electric4T.Impp(:,2))./4.55; %raltive to STC

bplot1 = plot(my_ax1,mismatch);%edgeless bar plot




%Adjusting left axes
box(my_ax1,'on');
my_ax1.LineWidth = 1.1;
my_ax1.YLim = [1.1*min(mismatch) 1.1*max(mismatch)];
my_ax1.XLim = [0 8760];
my_ax1.XTick = monthsInHours;%only write the ticks where the bars are
my_ax1.XTickLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};%x-ticks labels
my_ax1.XTickLabelRotation = 45;%rotate the x-ticks
my_ax1.XLabel.String = 'Month';
my_ax1.XLabel.FontWeight = 'Bold';
my_ax1.YLabel.String = 'Current mismatch rel. to STC [%]';
my_ax1.YLabel.FontWeight = 'Bold';
my_ax1.FontName = 'Calibri';
my_ax1.FontSize = 14;



%Legends - Only show the markers
%lgd = legend(my_ax1,{'2 Terminal', '3 Terminal'}, 'Location','northeast');
%lgd.Box = 'off';%Hide the box around the legend

my_ax1.FontName = 'Calibri';
my_ax1.FontSize = 14;

%Saving the image
print(my_fig1,'CurrenmisMatch','-dmeta')