function S = SettingsLCOE_PVplant

%Characteristics battery
S.eff_round_trip=0.98;
S.SOC_initial=0.95;
S.SOC_max=0.95;
S.SOC_min=0.05;

%Characteristics for prices
S.Ins=0.02; S.M=0.01;S.n=20;S.g=0.03;S.e=0.02;S.e_eng=0.03;S.d=0.03; %NPC details for economic analysis
S.hpf=3000; %€ fixed cost for heat pump 
S.hpv=500; %€/kW variable cost for heat pump
S.pvf=2000; %€ fixed cost for pv
S.pvv=1500; %€/kW variable cost for pv
S.bf=600; %€ fixed cost for battery
S.bv=550; %€/kWh variable cost for battery
S.pv_size=1.08;

% Average and amplitude for bought electricity cost
S.average_bought_cost = 0.25;  % €/kWh
S.amplitude_bought_cost = 0.02;  % €/kWh

% Average and amplitude for sold electricity cost
S.average_sold_cost = 0.06;  % €/kWh
S.amplitude_sold_cost = 0.005;  % €/kWh

end