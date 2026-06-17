function [LCOE_output, TOOLBOX_input] = LCOE_residentialPV(TOOLBOX_input,ELECTRIC_output,DEGRADATION_output,CONVERSION_output)
% LCOE_residential calculates the costs of the PV system in a residential building.
%
% This function calculates the NPC of a the PV system, including a heat
% pump and a battery storage
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% ELECTRIC_output : struct
%   Output of the ELECTRIC block
%
% Returns
% -------
% LCOE_output : struct
%   Output of this module block
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by Tekin Kodzhabash, integrated by Youri Blom


Pac = CONVERSION_output.Pac;
hours = 1:length(Pac);
if ~TOOLBOX_input.script
    TOOLBOX_input = get_LCOE_residential_input(TOOLBOX_input);
end

Size_HP = TOOLBOX_input.FinancialPart.Size_HP;
Size_Bat = TOOLBOX_input.FinancialPart.Size_Bat;

S = SettingsLCOE_PVplant;
%Characteristics battery
eff_round_trip=S.eff_round_trip;
eff_bat=sqrt(eff_round_trip);
SOC_initial=S.SOC_initial;
SOC_max=S.SOC_max;
SOC_min=S.SOC_min;

%Characteristics for prices
Ins=S.Ins; M=S.M;n=S.n; g=S.g;e=S.e;e_eng=S.e_eng;d=S.d; %NPC details for economic analysis
hpf=S.hpf; %€ fixed cost for heat pump 
hpv=S.hpv; %€/kW variable cost for heat pump
pvf=S.pvf; %€ fixed cost for pv
pvv=S.pvv; %€/kW variable cost for pv
bf=S.bf; %€ fixed cost for battery
bv=S.bv; %€/kWh variable cost for battery
pv_size=S.pv_size;

% Average and amplitude for bought electricity cost
average_bought_cost = S.average_bought_cost;  % €/kWh
amplitude_bought_cost = S.amplitude_bought_cost;  % €/kWh
% Equation for the cost of bought electricity
bought_cost = average_bought_cost + amplitude_bought_cost * sin((pi / 12) * (hours - 6));

% Average and amplitude for sold electricity cost
average_sold_cost = S.average_sold_cost;  % €/kWh
amplitude_sold_cost = S.amplitude_sold_cost;  % €/kWh

% Equation for the cost of sold electricity
sold_cost = average_sold_cost + amplitude_sold_cost * sin(pi/6*(hours-5));



bat_max=SOC_max*Size_Bat;
bat_min=SOC_min*Size_Bat;

[electricity_consumption, boilerWall] = heatpumpint(TOOLBOX_input,Size_HP);

bat_state=zeros(length(hours)+1,1);
bat_state(1) = SOC_initial*Size_Bat;
grid=zeros(length(hours),1);
E_bought=zeros(length(hours),1);
E_sold=zeros(length(hours),1);
for t=hours
    grid(t)=electricity_consumption(t)-Pac(t);
    if Pac(t)>=electricity_consumption(t) %There is still available ac power after the electricity consumption of the house, INVERTER SHOULD BE UPDATED
        %Charge the battery with excess PV or cell to the grid
        if bat_state(t)<bat_max %first check if the battery needs to be charged
            bat_state(t+1)=min(bat_max,-grid(t)*eff_bat+bat_state(t)); % Battery can only be charged up to maximum limit
            if (bat_state(t+1)-bat_state(t))/eff_bat+electricity_consumption(t)<=Pac(t) %Battery is charged but there is still excess energy
                E_sold(t)=Pac(t)-(((bat_state(t+1)-bat_state(t))/eff_bat)+electricity_consumption(t));
            else %battery is charged and there is no excess energy to be sold
                E_sold(t)=0;
            end

        else %Battery is already at the maximum capacity, electricity is sold to the grid
            E_sold(t)=-grid(t);
            bat_state(t+1)=bat_max;

        end
    else %Ac power of PV is not enough to cover the electricity consumption, or there is no production at all
        if bat_state(t)>bat_min %Battery can discharge to cover the electricity consumption of the house %SOMETHING NEEDS TO BE DONE HERE
            bat_state(t+1)=max(bat_min,bat_state(t)-grid(t)./(eff_bat)); % Battery can be discharged up to its minimum limit
            if bat_state(t)-grid(t)/(eff_bat)<bat_min %Discharging of battery is not enough to cover the electricity consumption of the house
                E_bought(t)=-(bat_state(t)-bat_state(t+1)).*eff_bat+grid(t);

            else
                E_bought(t)=0;
            end

        else %Battey at the minumum state, therefore electricity is bought from the grid
            E_bought(t)=grid(t);%Ac power is inclueded to consider the self consumption of the inverter
            bat_state(t+1)=bat_min;
        end
    end


end

% Net present cost calculation
bat_cost=bf+bv.*Size_Bat; if Size_Bat ==0; bat_cost = 0; end
heatpump_cost=hpf+hpv.*Size_HP; if Size_HP ==0; heatpump_cost = 0; end
I_o=bat_cost + heatpump_cost+pvf+pvv*pv_size;
Ic=pvf+pvv*pv_size+heatpump_cost; %will change with pv size

C_sold = sum(E_sold.*sold_cost');
C_bought = sum(E_bought.*bought_cost');

bat_cost_replacement = (bat_cost ./ ((1 + d)^15)-(bat_cost.*2)./(3*(1+d)^20))';

C_boiler=boilerWall*0.1253/(1000*0.92);
C_boiler_sum=sum(C_boiler);

syms k
SUM_op=symsum(((C_bought-C_sold).*(1+e).^k+C_boiler_sum.*(1+e_eng).^k+(M+Ins).*Ic.*(1+g).^k)/(1+d)^k,k,1,20);
SUM_op_result=eval(SUM_op);

NPC = I_o + bat_cost_replacement+ SUM_op_result;

LCOE_output.NPC = NPC;
LCOE_output.E_bought = E_bought;
LCOE_output.E_sold = E_sold;
LCOE_output.bat_state = bat_state;
LCOE_output.grid = grid;

disp(append('The NPC of this system is: ',num2str(round(NPC))))
end