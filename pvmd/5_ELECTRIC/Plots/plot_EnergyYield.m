function plot_EnergyYield(A, Pmpp, time,period,TOOLBOX_input)
%plot_EnergyYield plots the energy yield over the different months
%
% This function makes a figure of the produced electricity over the whole
% year
%
% Parameters
% ----------
% A : double
%   Absorbed irradiance
% Pmpp : double
%   The maximum powerpoint at each hour
% time : double
%   Time of which the simulation is performed
% period : double
%   Period of which the simulation is performed
% TOOLBOX_input : struct
%   Inputs of the simulation
%
% Returns
% -------
%(none)
%
% Developed by by Abdallah Nour El Din and Malte Vogt


if TOOLBOX_input.runPeriodic
    if time <= 168
        figure(71)
        hold on
        for i=1:time
            i1=area([i-0.5 i+0.5],...
                [mean(A(i,1:end)) mean(A(i,1:end))],...
                'Linestyle','none','FaceColor','y');
            p1=area([i-0.5 i+0.5],...
                [Pmpp(i) Pmpp(i)],...
                'Linestyle','none','FaceColor','b');
            alpha(0.66);
        end
        hold off
        grid on
        ylabel('Energy [Wh]','FontSize',14);
        xlabel('Time [hours]','FontSize',14);
        legend([i1 p1],'Incident Irradiance','DC Power Output');
    else
        figure(72)
%         clo(figure(72))
        hold on
        Pdaily=zeros(time/24,1);
        Months=[31,28,31,30,31,30,31,31,30,31,30,31];
        cell1=sum(Months(1:period(2)-1))*24+(period(1)-1)*24;
        for i=1:time/24
            Pdaily(i)=sum(Pmpp((i-1)*24+1:i*24));
            area([cell1/24+i-0.5 cell1/24+i+0.5],...
                [mean(sum(A((i-1)*24+1:i*24,1:end))), ...
                mean(sum(A((i-1)*24+1:i*24,1:end)))],'Linestyle','none','FaceColor','y');
            area([cell1/24+i-0.5 cell1/24+i+0.5],...
                [Pdaily(i) Pdaily(i)],...
                'Linestyle','none','FaceColor','b');
            alpha(0.66);    
        end
        hold off
        grid on
        ylabel('DC Energy Output [Wh]','FontSize',14);
        if time/24 >=31
            xticks([15.5 31 45 59 74.5 90 ...
                105 120 135.5 151 166 181 197 ...
                212 227.5 243 258 273 288 303 317 332 350])
            xticklabels({'Jan','',...
                'Feb','','Mar','','Apr','','May','','Jun',...
                '','Jul','','Aug','','Sep','','Oct',...
                '','Nov','','Dec'})
        else
            xlabel('Time [Days]')
        end
        grid on
    ylabel('Daily Energy Yield [Wh]')
    legend('Incident Power','DC Power Output');
    title('Incident Energy vs System DC Power Output');
    axis tight;
    grid on;
    end
% xlim([0 time/24])
else %non-periodic simulations
    if time <= 168
        figure(71)
        for j = 1:numel(A)
            subplot(size(A,1),size(A,2),j);
            for i=1:time
                area([i-0.5 i+0.5],...
                    [mean(A{j}(i,1:end)) mean(A{j}(i,1:end))],...
                    'Linestyle','none','FaceColor','y');
                hold on;
                area([i-0.5 i+0.5],...
                    [Pmpp{j}(i) Pmpp{j}(i)],...
                    'Linestyle','none','FaceColor','b');
                alpha(0.66);
            end
            grid on
            ylabel('Energy [Wh]','FontSize',14);
            xlabel('Time [hours]','FontSize',14);
            legend(gca,{'Incident Irradiance','DC Power Output'});
            title("DC Output & Irradiance Plot, Module: " + j);
        end  
    else
        figure(72)
%         clo(figure(72))
        Pdaily=zeros(time/24,1);
        Months=[31,28,31,30,31,30,31,31,30,31,30,31];
        cell1=sum(Months(1:period(2)-1))*24+(period(1)-1)*24;
        for j = 1:numel(A)
            subplot(size(A,1),size(A,2),j);
            hold on
            for i=1:time/24
                Pdaily(i)=sum(Pmpp{j}((i-1)*24+1:i*24));
                area([cell1/24+i-0.5 cell1/24+i+0.5],...
                    [mean(sum(A{j}((i-1)*24+1:i*24,1:end))), ...
                    mean(sum(A{j}((i-1)*24+1:i*24,1:end)))],'Linestyle','none','FaceColor','y');
                area([cell1/24+i-0.5 cell1/24+i+0.5],...
                    [Pdaily(i) Pdaily(i)],...
                    'Linestyle','none','FaceColor','b');
                alpha(0.66);    
            end
            hold off
            grid on
            ylabel('DC Energy Output [Wh]','FontSize',14);
            if time/24 >31
                xticks([15.5 31 45 59 74.5 90 ...
                    105 120 135.5 151 166 181 197 ...
                    212 227.5 243 258 273 288 303 317 332 350])
                xticklabels({'Jan','',...
                    'Feb','','Mar','','Apr','','May','','Jun',...
                    '','Jul','','Aug','','Sep','','Oct',...
                    '','Nov','','Dec'})
            else
                xlabel('Time [Days]')
            end
            title("DC Power Output, Module: " + j);
        end
    end    
end
