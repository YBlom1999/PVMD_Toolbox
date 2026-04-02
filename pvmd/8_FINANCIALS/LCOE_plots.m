function LCOE_plots(scenario,Cost_1,Cost_2,Cost_3)

switch nargin
    case 2
        x1=1:length(Cost_1);
        plot(x1,Cost_1,'r')
        legend(sprintf('Scenario %d',scenario))
        title('LCOE(EUR per kWp) vs Expected time period')
        xlabel('Time period(in years)')
        ylabel('LCOE(EUR per kWp) values')
        grid minor;
    case 4
        x1=1:length(Cost_1);
        x2=1:length(Cost_2);
        x3=1:length(Cost_3);
        hold on
        plot(x1,Cost_1,'k')
        plot(x2,Cost_2,'r')
        plot(x3,Cost_3,'g') 
        hold off 
        legend('Same module till expected period','Same technology reinstalled at Lifetime','New technology installed at Lifetime')
        title('LCOE(EUR per kWp) vs Expected time period')
        xlabel('Time period(in years)')
        ylabel('LCOE(EUR per kWp) values')
   