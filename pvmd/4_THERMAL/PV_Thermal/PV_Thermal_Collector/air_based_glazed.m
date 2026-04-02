function [pvt_collector_output] = air_based_glazed(MODULE_output, WEATHER_output)
% Calculates the thermal and PV yield of glazed air-based PVT collector
% 06-23
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
w = W*0.9375;
Ia(Ia<=0.05) = 0.05;

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

pvt_collector_output.y = zeros(length(Ia),6);

x0 = [25, 24, 25, 25, 24, 23];
h = 24/length(Ia):24/length(Ia):24;
tspan = 0:1:length(Ta)-1;
dm_f = 0.02;                     % mass flow rate (kg/sec)
[t1, x1] = ode23(@(t, x) odeq(t, x, Ta, G, V_w, L, W, w, dm_f), tspan, x0);
  
T_amb = Ta(floor(t1) + 1);
T_fin = T_amb;

d = 0.02;
A_mod = L*W;                          % Module Area 1.28x.32 m2
D_h = 2*(w*d)/(w+d);                  % Hydraulic Diameter 
lamda_f = 0.026;                      % Thermal conductivity of fluid
C_p = 1007;                           % Heat capacity (J/kgK)
mu_f = 1.849e-5;
rou = 1.184;
V = 1.75;
I_ref = 1000;
Tcel_ref = 25;
beta_p = 0.0045;                      % Temperature coefficient
delta = 0.052;                        % Solar radiation coefficient
eta_ref = 0.15;

Re = (rou*V*D_h/mu_f);                % Reynolds number
Pr = (mu_f*C_p)/lamda_f;              % Prandtl number
Nu = 0.023*(Re^0.8)*(Pr^0.4);         % Nusselt number

hvfmt = Nu*((lamda_f)/D_h);           % Convective
hvfmi = Nu*((lamda_f)/D_h);           % Convective

%% for air : Tf
E1a = x1(:,5).'*hvfmt+ x1(:,6).'*hvfmi;
E2a = hvfmt + hvfmi;
EE = E1a/E2a;
EEE = - w*E2a/(dm_f*C_p);
Tfo_1 = (T_fin - EE)*exp(EEE*L) + EE;
Tmean = (T_fin + Tfo_1)/2;
% y = x1;               % to see results on command window

Ttg_1 = x1(:, 1)';      % Store Ttg values in the first column of y
Tcf_1 = x1(:, 2)';      % Store Tcf values in the second column of y
Tbg_1 = x1(:, 3)';      % Store Tbg values in the third column of y
Tc_1 = x1(:, 4)';       % Store Tc values in the fourth column of y
Tt_1 = x1(:, 5)';       % Store Tt values in the fifth column of y
Ti_1 = x1(:, 6)';       % Store Ti values in the sixth column of y

groupSize = 3600; % Number of elements in each group
% Calculate the number of groups
numGroups = floor(numel(Tfo_1)/groupSize);
% Calculate the average of each column (group)
Two = mean(reshape(Tfo_1(1:numGroups*groupSize), groupSize, numGroups));
Tcf = mean(reshape(Tcf_1(1:numGroups*groupSize), groupSize, numGroups));
T_fin_h = mean(reshape(T_fin(1:numGroups*groupSize), groupSize, numGroups));
Ttg = mean(reshape(Ttg_1(1:numGroups*groupSize), groupSize, numGroups));
Tbg = mean(reshape(Tbg_1(1:numGroups*groupSize), groupSize, numGroups));
Tc = mean(reshape(Tc_1(1:numGroups*groupSize), groupSize, numGroups));
Tt = mean(reshape(Tt_1(1:numGroups*groupSize), groupSize, numGroups));
Ti = mean(reshape(Ti_1(1:numGroups*groupSize), groupSize, numGroups));

eta_E = eta_ref.*(1 - (beta_p.*(Tc - Tcel_ref))+(delta.*log(Ia/I_ref)));
eta_E(Ia<=1) = 0;
eta_tha = dm_f.*C_p.*(Two - T_fin_h)./(A_mod.*Ia);
eta_tha(eta_tha<0) = 0;
eta_tha(eta_tha>0.5) = 0.5;
eta_tha(Ia<=1) = 0;

Quele = Ia.*eta_E;
Quth = Ia.*eta_tha;

% Total solar irradiance received by the PVT system in kWh/m^2/year
total_I_sun = sum(mean(WEATHER_output.Irr,2))./1000; 
% Average thermal efficiency
etaSum = sum(eta_tha)./length(I1);
% Average electrical efficiency
etaE = sum(eta_E)/length(I1);

% Total thermal output in kWh/m2
total_thermal_output = etaSum*total_I_sun;
% Total electrical output in kWh/m2
total_electrical_output = etaE*total_I_sun;

pvt_collector_output.eta_tha=eta_tha;
pvt_collector_output.etaSum=etaSum;
pvt_collector_output.etaE=etaE;
pvt_collector_output.eta_E=eta_E;
pvt_collector_output.Two=Two;
pvt_collector_output.Tcf=Tcf;
pvt_collector_output.Ttg=Ttg;
pvt_collector_output.Tbg=Tbg;
pvt_collector_output.T=Tc;
pvt_collector_output.Tt=Tt;
pvt_collector_output.Ti=Ti;
pvt_collector_output.Tmean=Tmean;
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

f = figure(101);
plot(xx, T_fin_h,'DisplayName','T_a')
hold on
plot(xx, Ttg,'DisplayName','T_{tg}')
plot(xx, Tcf,'DisplayName','T_{cf}')
plot(xx, Tbg,'DisplayName','T_{bg}')
plot(xx, Tc,'DisplayName','T_{c}')
plot(xx, Tt,'DisplayName','T_{t}')
plot(xx, Ti,'DisplayName','T_{i}')
plot(xx, Two,'DisplayName','T_{f,o}')
ylabel('temperature (^oC)');
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
%     legend('$T_a$','$T_{\mathcal{G}}$','$T_{cf}$','$T_{\mathcal{C}}$','$T_{c}$','$T_{t}$','$T_{i}$','$T_{f,o}$','Interpreter','LaTeX','Fontsize',12)
    SetFigureDefaults(f)

    f = figure(1012);
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
    title('Power Output Plot');
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

function X = odeq(t, x, Ta, G, V_w, L, W, w, dm_f)
t2 = floor(t + 1);
Is = G(t2);
T_amb = Ta(t2);
T_fin = T_amb;
V_w2 = V_w(t2);

V = 1.75;
d = 0.01;
A_mod = L*W;                     % Module Area 1.28x.32 m2
A_c = L*w;                       % Horizontal Area of Air duct 1.28x.3 m2
alpha_gla = 0.1;                 % Absorptivity of glass cover
C_gla = 670;               
sigma = 5.67e-8;
em_gla = 0.93;                   % Emissivity of glass cover
l_gla = 0.0032;                  % Thickness of glass cover (m)
rou_gla = 2200;
M_gla = rou_gla*l_gla*A_mod;
l_cel = 0.0003;                  % Thickness of solar cell (m)
lamda_gla = 1.1;                 % thermal conductivity of glass cover (W/mK)
lamda_cel = 148;                 % thermal conductivity of solar cell (W/mK)
rou_cel = 2330;
M_cel = rou_cel*l_cel*A_c;
C_cel = 900;
tau_gla = 0.84;                  % Transmisivity of glass cover
alpha_cel = 0.90;                % Absorptivity of solar cell
beta = 0.88;                     % Packing factor
l_ted = 0.00075;                 % Thickness of Tedlar (m)
rou_ted = 1200;
M_ted = rou_ted*l_ted*A_mod;
C_ted = 1200;
alpha_ted = 0.8;                 % Absorptivity of tedlar
zeta = 0.984;                    % Report duct area/PV module area
em_ins = 0.6;                    % Emissivity of insulator
em_ted = 0.85;                   % Emissivity of tedlar
lamda_ted = 0.2;                 % thermal conductivity of Tedlar (W/mK)
D_h = 2*(w*d)/(w+d);             % Hydraulic Diameter 
mu_f = 1.5e-5;                   % Kinematic viscosity of fluid (m2/sec)
lamda_f = 0.026;
% lamda_ins = 0.034;             % thermal conductivity of insulator (W/mK)
C_p = 1007;
rou_ins = 20;
l_ins = 0.04;
M_ins = rou_ins*l_ins*A_mod;
C_ins = 670;
beta_p = 0.0045;                 % Temperature coefficient
delta = 0.052;                   % Solar radiation coefficient
eta_ref = 0.15;                  % Reference Efficiency
Tcel_ref = 25;
I_ref = 1000;
rou = 1.184;
delta_a = 0.02;
M_f = rou*A_mod*delta_a;
beta_t = 3.43e-3;                % coefficient of thermal expansion for air
g = 9.81;

Tsk = 0.0552* (T_amb+273.15)^1.5;
Tsky = Tsk - 273.15;             % Temperature of Sky
% hva = 6.5 + 3.3*V_w2;          % Convective heat coefficient due to wind
hva(V_w2<5) = 5.7 + 3.8*V_w2;
hva(V_w2>=5) = 6.47 + V_w2^0.78;

Re = (rou*V*D_h/mu_f);           % Reynolds number
Pr = (mu_f*C_p)/lamda_f;         % Prandtl number
Nu = 0.023*(Re^0.8)*(Pr^0.4);    % Nusselt number

hcgmc = 1/((l_gla/lamda_gla) + (l_cel/lamda_cel));
hccmg = hcgmc;
hccmt = 1/((l_cel/lamda_cel) + (l_ted/lamda_ted));
hvfmt = Nu*((lamda_f)/D_h);
hvfmi = Nu*((lamda_f)/D_h);

%% for air : Tf
E1a = x(5)*hvfmi + x(6)*hvfmt;
E2a = hvfmi + hvfmt;

EE = E1a/E2a;
EEE = - w*E2a/(dm_f*C_p);
Tfo = (T_fin - EE)*exp(EEE*L) + EE;
Tf = (Tfo + T_fin)/2;
%% for glass : x(1)
alpha_d = lamda_f/(rou*C_p);     % thermal diffusivity - m²/s
Diff_1 = x(1) - x(2);
Diff_1(Diff_1<=0) = 0.001;
R_a_tg = (g*beta_t*(delta_a^3)*(Diff_1))/(mu_f*alpha_d);     % Diff_1 = x(1) - x(2)
Diff_2 = x(3) - x(2);
Diff_2(Diff_2<=0) = 0.001;
R_a_bg = (g*beta_t*(delta_a^3)*(Diff_2))/(mu_f*alpha_d);     % Diff_2 = x(3) - x(2)

A_tgone = (1 - (1708)/(R_a_tg*cos(pi/6)));
A_tgone = max(A_tgone, 0);       % Take only positive values by setting negative values to zero
A_tgtwo = (((R_a_tg*cos((pi/6))/5830)^(1/3)) - 1);
A_tgtwo = max(A_tgtwo, 0);       % Take only positive values by setting negative values to zero
Nua_tg = (1 + 1.44*(A_tgone*(1 - (1708*(sin(1.8*(pi/6)))^1.6)/(R_a_tg*cos((pi/6))))) + A_tgtwo);

B_bgone = (1 - (1708)/(R_a_bg*cos(pi/6)));
B_bgone = max(B_bgone, 0);       % Take only positive values by setting negative values to zero
B_bgtwo = (((R_a_bg*cos((pi/6))/5830)^(1/3)) - 1);
B_bgtwo = max(B_bgtwo, 0);       % Take only positive values by setting negative values to zero
Nua_bg = (1 + 1.44*(B_bgone*(1 - (1708*(sin(1.8*(pi/6)))^1.6)/(R_a_bg*cos((pi/6))))) + B_bgtwo);

hvtgcf = (Nua_tg*lamda_f)/delta_a;
hvcfbg = (Nua_bg*lamda_f)/delta_a;
hrtgmbg = sigma*((x(1) + x(3))*(x(1)^2 + x(3)^2))/((1/em_gla)+(1/em_gla) - 1);
hrtgms = sigma*em_gla*(x(1)^2 + Tsky^2)*(x(1) + Tsky);

%% for top glass : x(1)
x1 = (alpha_gla*Is + hrtgms*(Tsky - x(1)) + hva*(T_amb - x(1)) - hrtgmbg*(x(1) - x(3)) - hvtgcf*(x(1) - x(2)))*(A_mod/(M_gla*C_gla));

%% for closed air space : x(2)
x2 = (hvcfbg*(x(3) - x(2)) - hvtgcf*(x(2) - x(1)))*(A_mod/(M_f*C_p));

%% for glass below : x(3)
x3 = ((alpha_gla^2)*Is - hrtgmbg*(x(3) - x(1)) + hvcfbg*(x(2) - x(3)) - hcgmc*(x(3) - x(4)))*(A_mod/(M_gla*C_gla));

%% for PV solar cell : x(4)
Quele = Is*A_mod*eta_ref*(1 - (beta_p*(x(4) - Tcel_ref))+(delta*log(Is/I_ref)));
x4 = ((tau_gla^2)*alpha_cel*Is*beta - hccmg*(x(4) - x(3)) - hccmt*(x(4) - x(5)))*(A_mod/(M_cel*C_cel)) - (Quele/(M_cel*C_cel));

%% for tedlar : x(5)
hrtmi = sigma*((x(6) + x(5))*(x(6)^2 + x(5)^2))/((1/em_ins)+(1/em_ted) - 1);
hrimt = hrtmi;
x5 = ((tau_gla^2)*alpha_ted*Is*(1 - beta) + hccmt*(x(4) - x(5)) - hvfmt*zeta*(x(5) - Tf) - hrimt*zeta*(x(5) - x(6)))*A_mod/(M_ted*C_ted);

%% for insulator : x(6)
x6 = (hrtmi*(x(5) - x(4)) + hvfmi*(Tf - x(6)) - hva*(x(6) - T_amb))*A_c/(M_ins*C_ins);

X  = [x1; x2; x3; x4; x5; x6];
end
