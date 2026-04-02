function [pvt_collector_output] = water(MODULE_output, WEATHER_output)
% Calculates the thermal and PV yield of water-based PVT collector
% 08-23 - Water-based PVT: ths
%
% Parameters
% ----------
% MODULE_output : struct
%   Output of this module block
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% pvt_collector_output : struct
%   Simulation results of the PVT module
% 
% Developed by ZUA.

Ia = mean(WEATHER_output.Irr,2)';
Tam = [WEATHER_output.ambient_temperature]';
Vw = [WEATHER_output.wind_speed]';
L = MODULE_output.ML;
W = MODULE_output.MW;
d_o = 0.015;
d_i = 0.013;            
dm_w = 0.01;                          % mass flow rate (kg/sec)
Ia(Ia<=0.05) = 0.05;                  % so that efficiency doesn't give NaN values

idx = find(Ia>0.05);                  % find the indices of non-zero elements in Ia
if length(idx) > length(Ia)           % if there are non-zero elements in Ia
    idx = idx(1:length(Ia));          % select the non-zero elements
end
I1 = Ia(idx);                         % length when there is solar radiation > 0.05  

% Define the number of times to repeat Ta and Ia, since dynami model so each hour is converted in sec
n_repeats = 3600;
% Repeat Ta and G n_repeats times
Ta = repelem(Tam, n_repeats);
G = repelem(Ia, n_repeats);
V_w = repelem(Vw, n_repeats);

pvt_collector_output.y = zeros(length(Ia),7);

x0 = [25, 25 25, 25, 25, 25, 20];
h = 24/length(Ia):24/length(Ia):24;
tspan = 0:1:length(Ta)-1;
[t1, x1] = ode23(@(t, x) odeq(t, x, Ta, G, V_w, L, W, d_o, d_i, dm_w), tspan, x0);

T_amb = Ta(floor(t1) + 1);
Twin = T_amb;

Two_1 = 2*x1(:,6).' - Twin;

A_mod = L*W;                     % Module Area
beta_p = 0.0045;                 
delta = 0.052;                   
eta_ref = 0.15;                 
Tcel_ref = 25;
I_ref = 1000;
C_w = 4186;

Tg_1 = x1(:,1).';
Tc_1 = x1(:,2).';
Tt_1 = x1(:,3).';
Tr_1 = x1(:,4).';
Tm_1 = x1(:,5).';
Tw_1 = x1(:,6).';
Ti_1 = x1(:,7).';

groupSize = 3600; % Number of elements in each group
% Calculate the number of groups
numGroups = floor(numel(Two_1)/groupSize);
% Calculate the average of each column (group)
Two = mean(reshape(Two_1(1:numGroups*groupSize), groupSize, numGroups));
Tg = mean(reshape(Tg_1(1:numGroups*groupSize), groupSize, numGroups));
Tc = mean(reshape(Tc_1(1:numGroups*groupSize), groupSize, numGroups));
Tt = mean(reshape(Tt_1(1:numGroups*groupSize), groupSize, numGroups));
Tr = mean(reshape(Tr_1(1:numGroups*groupSize), groupSize, numGroups));
Tm = mean(reshape(Tm_1(1:numGroups*groupSize), groupSize, numGroups));
Tw = mean(reshape(Tw_1(1:numGroups*groupSize), groupSize, numGroups));
Ti = mean(reshape(Ti_1(1:numGroups*groupSize), groupSize, numGroups));
Twin_h = mean(reshape(Twin(1:numGroups*groupSize), groupSize, numGroups));

eta_E = eta_ref.*(1 - (beta_p.*(Tc - Tcel_ref))+(delta*log(Ia/I_ref)));
eta_E(Ia<=1) = 0;
eta_tha = dm_w.*C_w.*(Two - Twin_h)./(A_mod.*Ia);
eta_tha(eta_tha<0) = 0;
eta_tha(eta_tha>0.5) = 0.5;
eta_tha(Ia<=1) = 0;

Quele = Ia.*eta_E;
Quth = Ia.*eta_tha;

pvt_collector_output.y(:, 1) = Tg;        % Store Tg values in the first column of y
pvt_collector_output.y(:, 2) = Tc;        % Store Tc values in the second column of y
pvt_collector_output.y(:, 3) = Tt;        % Store Tt values in the third column of y
pvt_collector_output.y(:, 4) = Tr;        % Store Tt values in the fourth column of y
pvt_collector_output.y(:, 5) = Tm;        % Store Tt values in the fifth column of y
pvt_collector_output.y(:, 6) = Tw;        % Store Two values in the sixth column of y
pvt_collector_output.y(:, 7) = Ti;        % Store Two values in the seventh column of y
pvt_collector_output.y(:, 8) = Two;       % Store Two values in the eight column of y
pvt_collector_output.y(:, 9) = eta_tha;   % Store eta_tha values in the ninth column of y
pvt_collector_output.y(:, 10) = eta_E;    % Store eta_E values in the tenth column of y

% Total solar irradiance received by the PVT system in kWh/m^2/year
total_I_sun = sum(mean(WEATHER_output.Irr,2))/1000; 
% Average thermal efficiency
etaSum = sum(eta_tha)/length(I1);
% Average electrical efficiency
etaE = sum(eta_E)/length(I1);

%% change if kWh or kWh/m2
% Yearly total thermal output in kWh/m2
total_thermal_output = etaSum*total_I_sun;
% Yearly total electrical output in kWh/m2
total_electrical_output = etaE*total_I_sun;

pvt_collector_output.eta_tha=eta_tha;
pvt_collector_output.etaSum=etaSum;
pvt_collector_output.etaE=etaE;
pvt_collector_output.eta_E=eta_E;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.total_thermal_output=total_thermal_output;
pvt_collector_output.total_electrical_output=total_electrical_output;
pvt_collector_output.h=h;

%% Plot
period=WEATHER_output.Period;
Irr = mean(WEATHER_output.Irr,2);
ambient_temperature = WEATHER_output.ambient_temperature;

time = length(Irr);
if time <= 24
    x = 1:time;
    % Irradiance
    f = figure(33); hold on;
 
    yyaxis right
    fill([x, fliplr(x)],[Irr; zeros(time,1)],'y','LineStyle','none')
    alpha(0.33)
    ylabel('Irradiance [W/m^2]');
    f.CurrentAxes.YColor = 'k';
    
    yyaxis left
    plot(x,ambient_temperature,'b-','LineWidth',2);
    f.CurrentAxes.YColor = 'k';
    hold off

    xlabel('Time [Hours]'); ylabel('Temperature [^oC]');
    legend({'T_{amb}','Irradiance'},'Location','NorthEast',...
        'FontSize',14);
    xlim([1 time])
    SetFigureDefaults(f)
else
    xx = 1:time;

% Xticks settings (per hour data)
    months = 24*[31,28,31,30,31,30,31,31,30,31,30,31];  
    ticks = 24*[15.5 31 45 59 74.5 90 105 120 135.5 151 166 181 197 212 227.5 ...
        243 258 273 288 303 317 332 350]-(period(1)+sum(months(1:period(2)-1)));
    month_labelss = {'Jan','','Feb','','Mar','','Apr','','May','','Jun',...
            '','Jul','','Aug','','Sep','','Oct','','Nov','','Dec'};
end

f = figure(100);
plot(xx,Tam,'DisplayName','T_{a}')
hold on
plot(xx, Tg,'DisplayName','T_{gla}')
plot(xx, Tc,'DisplayName','T_{cel}')
plot(xx, Tt,'DisplayName','T_{ted}')
plot(xx, Tr,'DisplayName','T_{abs}')
plot(xx, Tm,'DisplayName','T_{mt}')
plot(xx, Two,'DisplayName','T_{w,o}')
plot(xx, Tw,'DisplayName','T_{w}')
plot(xx, Ti,'DisplayName','T_{ins}')
ylabel('Temperature [^oC]');
yyaxis right
plot(xx, Ia,'DisplayName','I_{sun}')
ylabel('Solar radiation [Wm^{-2}]');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Temperature Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

f = figure(1011);
plot(xx,Quele,'DisplayName','E_{p}')
ylabel('Electrical Power [Wm^{-2}]');
yyaxis right
plot(xx, Quth,'DisplayName','T_{p}')
ylabel('Thermal power [Wm^{-2}]');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Temperature Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

%%
function SetFigureDefaults(f)
    % Set several figure characteristics for correct representation
    % figure object to be modified

    f.CurrentAxes.FontSize = 14;
    f.CurrentAxes.Box = 'on';
    f.CurrentAxes.XGrid = 'on';
    f.CurrentAxes.YGrid = 'on';
    f.CurrentAxes.XLabel.FontSize = 16;
    f.CurrentAxes.YLabel.FontSize = 16;
    f.CurrentAxes.Title.FontSize = 16;
end

end

function X = odeq(t, x, Ta, G, V_w, L, W, d_o, d_i, dm_w)
t2 = floor(t + 1);
Is = G(t2);
T_amb = Ta(t2);
V_w2 = V_w(t2);
Twin = T_amb;

A_mod = L*W;                     % Module Area
A_c = A_mod*0.80;
alpha_gla = 0.1;                 % Absorptivity of glass cover
sigma = 5.67e-8;
l_gla = 0.0032;     
em_gla = 0.88;          
l_cel = 0.0003;           
lamda_gla = 1.1;          
lamda_cel = 148;          
tau_gla = 0.84;   % Slimani           
alpha_cel = 0.90;         
beta = 0.88;          
eta_ref = 0.15;         
beta_p = 0.0045;       
lamda_ted = 0.2;          
lamda_abs = 160;                 % for Aluminuml_i = 0.08;            
lamda_ins = 0.034;        
l_mo = 0.001;          
lamda_mo = 40;          
lamda_w = 0.6;          
rou_gla = 2200;
C_gla = 670;               
M_gla = rou_gla*l_gla*A_mod;       
rou_cel = 2330;
M_cel = rou_cel*l_cel*A_c;    
C_c = 900;             
l_ted = 0.00025;                 % Thickness of Tedlar (m)
rou_ted = 1200;         
M_ted = rou_ted*A_mod*l_ted;     
C_ted = 1200;  
rou_m = 8900;          
C_m = 385;                                  
rou_w = 998;                               
C_w = 4186;                                 
rou_ins = 20;                                
C_ins = 670;                                 
Tcel_ref = 25;
I_ref = 1000;
delta = 0.052;     
l_ins = 0.05;

A_tm = 4*d_o*L/2;
A_tr = (W - 4*d_o)*L;
A_rm = (4*d_o*L)/4;
A_ri = (W - 4*d_o)*L;
A_mi = 4*(d_o*L + (3.14*d_o*L/2));
A_w = 3.14*4*d_i*L;      

V_m = pi*L*((d_o^2) - (d_i^2));         
V_w = pi*L*(d_i^2);                   
V_i = A_mod*l_ins - (pi*L*(d_o^2));       
l_abs = 0.001;
rho_abs = 2710;
C_abs = 887;
M_abs = rho_abs*A_c*l_abs;
M_m = rou_m*V_m;                         
M_w = rou_w*V_w;                         
M_i = rou_ins*V_i;                       

hcgmc = 1/((l_gla/lamda_gla) + (l_cel/lamda_cel));
hccmt = 1/((l_cel/lamda_cel) + (l_ted/lamda_ted));
hctmm = lamda_ted/l_ted;
hctmr = 1/((l_ted/lamda_ted) + (l_abs/lamda_abs));
hcrmm = 2*lamda_abs/((W - d_o)/4);
hcrmi = lamda_ins/l_ins;
hcmmi = 2*lamda_ins/l_ins;
h_w = 4.364*lamda_w/d_i;
A_mw_hvmmw = 1/((1/(h_w*A_w)) + (l_mo/(lamda_mo*d_i*L)));

Tsky = T_amb-20;       
hva = 6.5 + 3.3*V_w2;           

%% for glass : x(1)
hrgms = sigma*em_gla*(x(1) + Tsky)*(x(2)^2 + Tsky^2);
x1 = (alpha_gla*Is - hrgms*(x(1) - Tsky) - hcgmc*(x(1) - x(2)) - hva*(x(1) - T_amb))*(A_mod/(M_gla*C_gla));

%% for PV solar cell : x(2)
Quele = Is*A_mod*eta_ref*(1 - (beta_p*(x(2) - Tcel_ref))+(delta*log(Is/I_ref)));
x2 = (tau_gla*alpha_cel*Is*beta + hcgmc*(x(1) - x(2)) - hccmt*(x(2) - x(3)))*(A_mod/(M_cel*C_c)) - (Quele/(M_cel*C_c));

%% for tedlar : x(3)
x3 = (A_mod*hccmt*(x(2) - x(3)) - A_tm*hctmm*(x(3) - x(5)) - A_tr*hctmr*(x(3) - x(4)))*(1/(M_ted*C_ted));

%% for absorber : x(4)
x4 = (A_tr*hctmr*(x(3) - x(4)) - A_rm*hcrmm*(x(4) - x(5)) - A_ri*hcrmi*(x(4) - x(7)))*(1/(M_abs*C_abs));

%% for tube : x(5)
x5 = (A_tm*hctmm*(x(3) - x(5)) + A_rm*hcrmm*(x(4) - x(5)) - A_mi*hcmmi*(x(5) - x(7)) - A_mw_hvmmw*(x(5) - x(6)))*(1/(M_m*C_m));

%% for outlet water : T_wo
Two = 2*x(6) - Twin;

%% for water : x(6)
x6 = (A_mw_hvmmw*(x(5) - x(6)) + dm_w*C_w*(Twin - Two))*(1/(M_w*C_w));

%% for insulator : x(7)
x7 = (A_mi*hcmmi*(x(5) - x(7)) + A_ri*hcrmi*(x(4) - x(7)) - A_mod*hva*(x(7) - T_amb))*(1/(M_i*C_ins));

%% state outputs
X  = [x1; x2; x3; x4; x5; x6; x7];

end
