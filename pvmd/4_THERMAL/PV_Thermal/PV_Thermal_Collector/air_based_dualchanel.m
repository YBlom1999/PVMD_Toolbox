function [pvt_collector_output] = air_based_dualchanel(MODULE_output, WEATHER_output)
% Calculates the thermal and PV yield of dual chanel air-based PVT collector
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

x0 = [23, 24, 25, 24, 25, 23];
h = 24/length(Ia):24/length(Ia):24;
tspan = 0:1:length(Ta)-1;
dm_f = 0.02;                     % mass flow rate (kg/sec)
[t1, x1] = ode23(@(t, x) odeq(t, x, Ta, G, V_w, L, W, w, dm_f), tspan, x0);

T_amb = Ta(floor(t1) + 1);
T_fin = T_amb;                        % for air-based inlet temp is taken same as ambient temp

d = 0.02;                             % depth of air duct
A_mod = L*W;                          % Module Area 1.28x.32 m2
D_h = 2*(w*d)/(w+d);                  % Hydraulic Diameter 
lamda_f = 0.026;                      % Thermal conductivity of fluid
C_p = 1007;                           % Heat capacity (J/kgK)
mu_f = 1.849e-5;                      % viscosity of air (kg/ms)
rou = 1.184;                          % density of air (kg/m3)
V = 1.75;
eta_ref = 0.15;                       % Reference Efficiency
beta_p = 0.0045;                      % Temperature coefficient
Tcel_ref = 25;                        % reference cell temp
I_ref = 1000;                         % reference solar radiation
delta = 0.052;                        % Solar radiation coefficient

Re = (rou*V*D_h/mu_f);                % Reynolds number
Pr = (mu_f*C_p)/lamda_f;              % Prandtl number
Nu = 0.023*(Re^0.8)*(Pr^0.4);         % Nusselt number

hvgaf = Nu*((lamda_f)/D_h);           % convection
hvggaf = Nu*((lamda_f)/D_h);          % convection
hvlfmt = Nu*((lamda_f)/D_h);          % convection
hvlfmr = Nu*((lamda_f)/D_h);          % convection
%% for lower air : Tlf
E1a = x1(:,4).'*hvlfmt + x1(:,5).'*hvlfmr;
E2a = hvlfmt + hvlfmr;
EE = E1a/E2a;
EEE = - w*E2a/(dm_f*C_p);
Tlfo_1 = (T_fin - EE)*exp(EEE*L) + EE;

%% for air above : Taf
E3a = x1(:,1).'*hvgaf + x1(:,2).'*hvggaf;
E4a = hvggaf + hvgaf;
EE_ = E3a/E4a;
EEE_ = - W*E4a/(dm_f*C_p);
Tafo_1 = (Tlfo_1 - EE_)*exp(EEE_*L) + EE_;
Tmean = (T_fin + Tafo_1)/2;
% y = x1;               % to see results on command window

Ttg_1 = x1(:, 1)';      % Store Ttg values in the first column of y
Tbg_1 = x1(:, 2)';      % Store Tbg values in the second column of y
Tc_1 = x1(:, 3)';       % Store Tc values in the third column of y
Tt_1 = x1(:, 4)';       % Store Tt values in the fourth column of y
Tr_1 = x1(:, 5)';       % Store Tr values in the fifth column of y
Ti_1 = x1(:, 6)';       % Store Ti values in the sixth column of y

groupSize = 3600; % Number of elements in each group
% Calculate the number of groups
numGroups = floor(numel(Tafo_1)/groupSize);
% Calculate the average of each column (group)
Tlfo = mean(reshape(Tlfo_1(1:numGroups*groupSize), groupSize, numGroups));
Two = mean(reshape(Tafo_1(1:numGroups*groupSize), groupSize, numGroups));
T_fin_h = mean(reshape(T_fin(1:numGroups*groupSize), groupSize, numGroups));
Ttg = mean(reshape(Ttg_1(1:numGroups*groupSize), groupSize, numGroups));
Tbg = mean(reshape(Tbg_1(1:numGroups*groupSize), groupSize, numGroups));
Tc = mean(reshape(Tc_1(1:numGroups*groupSize), groupSize, numGroups));
Tt = mean(reshape(Tt_1(1:numGroups*groupSize), groupSize, numGroups));
Tr = mean(reshape(Tr_1(1:numGroups*groupSize), groupSize, numGroups));
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
pvt_collector_output.Tlfo=Tlfo;
pvt_collector_output.Ttg=Ttg;
pvt_collector_output.Tbg=Tbg;
pvt_collector_output.T=Tc;
pvt_collector_output.Tt=Tt;
pvt_collector_output.Tr=Tr;
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

f = figure(102);
plot(xx, T_fin_h,'DisplayName','T_a')
hold on
plot(xx, Ttg,'DisplayName','T_{tg}')
plot(xx, Tbg,'DisplayName','T_{bg}')
plot(xx, Tc,'DisplayName','T_c')
plot(xx, Tt,'DisplayName','T_t')
plot(xx, Tr,'DisplayName','T_r')
plot(xx, Ti,'DisplayName','T_i')
plot(xx, Tlfo,'DisplayName','T_{lfo}')
plot(xx, Two,'DisplayName','T_{afo}')
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
    legend()
    SetFigureDefaults(f)

    f = figure(1013);
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
T_fin = T_amb;                   % for air-based inlet temp is taken same as ambient temp
V_w2 = V_w(t2);

V = 1.75;
d = 0.02;                        % depth of air duct
A_mod = L*W;                     % Module Area 1.28x.32 m2
A_c = L*w;                       % Horizontal Area of Air duct 1.28x.3 m2
alpha_gla = 0.1;                 % Absorptivity of glass cover
l_gla = 0.0032;                  % Thickness of glass cover (m)
C_gla = 670;                     % heat capacity of glass 
rou_gla = 2200;                  % density of glass
M_gla = rou_gla*A_mod*l_gla;     % mass
sigma = 5.67e-8;                 % Stefan–Boltzmann constant
em_gla = 0.93;                   % Emissivity of glass cover
l_cel = 0.0003;                  % Thickness of solar cell (m)
rou_cel = 2330;                  % cell density
lamda_gla = 1.1;                 % thermal conductivity of glass cover (W/mK)
lamda_cel = 148;                 % thermal conductivity of solar cell (W/mK)
C_cel = 900;
M_cel = rou_cel*A_c*l_cel;
tau_gla = 0.84;                  % Transmisivity of glass cover
alpha_cel = 0.90;                % Absorptivity of solar cell
beta = 0.88;                     % Packing factor
C_ted = 1200;
rou_ted = 1200;
l_ted = 0.00075;                 % Thickness of Tedlar (m)
M_ted = rou_ted*A_mod*l_ted;
alpha_ted = 0.8;                 % Absorptivity of tedlar
zeta = 0.984;                    % Report duct area/PV module area
% em_ins = 0.6;                  % Emissivity of insulator
em_ted = 0.85;                   % Emissivity of tedlar
lamda_ted = 0.2;                 % thermal conductivity of Tedlar (W/mK)
D_h = 2*(w*d)/(w+d);             % Hydraulic Diameter 
mu_f = 1.849e-5;                 % Dynamic viscosity of fluid (kg/m.sec)
lamda_f = 0.026;
lamda_ins = 0.034;               % thermal conductivity of insulator (W/mK)
C_p = 1007;
l_ins = 0.4;
C_ins = 670;
rou_ins = 20;
M_ins = rou_ins*A_mod*l_ins;
beta_p = 0.0045;                 % Temperature coefficient
delta = 0.052;                   % Solar radiation coefficient
eta_ref = 0.15;                  % Reference Efficiency
Tcel_ref = 25;
I_ref = 1000;
rou = 1.184;
% alpha_abs = 0.95;
em_abs = 0.25;
lamda_abs = 250;
l_abs = 0.002;
rho_abs = 2710;
C_abs = 887;
M_abs = rho_abs*A_c*l_abs;

Tsk = 0.0552* (T_amb+273.15)^1.5;
Tsky = Tsk - 273.15;             % Temperature of Sky
% hva = 6.5 + 3.3*V_w2;          % Convective heat coefficient due to wind
hva(V_w2<5) = 5.7 + 3.8*V_w2;
hva(V_w2>=5) = 6.47 + V_w2^0.78;

Re = (rou*V*D_h/mu_f);           % Reynolds number
Pr = (mu_f*C_p)/lamda_f;         % Prandtl number
Nu = 0.023*(Re^0.8)*(Pr^0.4);    % Nusselt number

hcggmc = 1/((l_gla/lamda_gla) + (l_cel/lamda_cel));
hccmgg = hcggmc;
hcrmi =  1/((l_abs/lamda_abs) + (l_ins/lamda_ins));
hcimr = hcrmi;
hccmt = 1/((l_cel/lamda_cel) + (l_ted/lamda_ted));
hvgaf = Nu*((lamda_f)/D_h);
hvggaf = Nu*((lamda_f)/D_h);
hvlfmt = Nu*((lamda_f)/D_h);
hvlfmr = Nu*((lamda_f)/D_h);
%% for lower air : Tlf
E1a = x(4)*hvlfmt + x(5)*hvlfmr;
E2a = hvlfmt + hvlfmr;

EE = E1a/E2a;
EEE = - w*E2a/(dm_f*C_p);
Tlfo = (T_fin - EE)*exp(EEE*L) + EE;
Tlf = (Tlfo + T_fin)/2;

%% for air above : Taf
E3a = x(1)*hvgaf + x(2)*hvggaf;
E4a = hvggaf + hvgaf;

EE_ = E3a/E4a;
EEE_ = - W*E4a/(dm_f*C_p);
Tafo = (Tlfo - EE_)*exp(EEE_*L) + EE_;
Taf = (Tlfo + Tafo)/2;

%% for top-glass : x(1)
hrgms = sigma*em_gla*(x(1)^2 + Tsky^2)*(x(1) + Tsky);
hrggg = sigma*((x(1) + x(2))*(x(1)^2 + x(2)^2))/((1/em_gla)+(1/em_gla) - 1);
x1 = ((alpha_gla*Is) + (hrgms*(Tsky - x(1))) + hrggg*(x(2) - x(1)) + (hva*(T_amb - x(1))) - (hvgaf*(x(1) - Taf)))*(A_mod/(M_gla*C_gla));

%% for bottom-glass : x(2)
x2 = ((alpha_gla^2*Is) + (hrggg*(x(1) - x(2))) + hvggaf*(Taf - x(2)) - (hcggmc*(x(2) - x(3))))*(A_mod/(M_gla*C_gla));

%% for PV solar cell : x(3)
Quele = Is*A_mod*eta_ref*(1 - (beta_p*(x(3) - Tcel_ref))+(delta*log(Is/I_ref)));
x3 = ((tau_gla^2*alpha_cel*Is*beta) - (hccmgg*(x(3) - x(2)) )-(hccmt*(x(3) - x(4))))*(A_mod/(M_cel*C_cel)) - (Quele/(M_cel*C_cel));

%% for tedlar : x(4)
hrrmt = sigma*((x(5) + x(4))*(x(5)^2 + x(4)^2))/((1/em_abs)+(1/em_ted) - 1);
hrtmr = hrrmt;
x4 = ((tau_gla^2*alpha_ted*Is*(1 - beta)) + (hccmt*(x(3) - x(4))) - (hvlfmt*zeta*(x(4) - Tlf)) - (hrrmt*zeta*(x(4) - x(5))))*A_mod/(M_ted*C_ted);

%% for absorber : x(5)
x5 = ((hcrmi*(x(6) - x(5))) - (hvlfmr*(x(5) - Tlf)) - (hrtmr*(x(5) - x(4))))*A_c/(M_abs*C_abs);

%% for insulator : x(6)
x6 = ((hcimr*(x(5) - x(6))) - (hva*(x(6) - T_amb)))*A_c/(M_ins*C_ins);

X  = [x1; x2; x3; x4; x5; x6];
end
