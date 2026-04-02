function [pvt_collector_output] = bifluid_air_water(MODULE_output, WEATHER_output)
% Calculates the thermal and PV yield of bi-fluid PVT collector
% 07-23
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

pvt_collector_output.y = zeros(length(Ia),9);
x0 = [23, 25, 24, 24, 23, 24, 23, 24, 20];
tspan = 0:1:length(Ta)-1;
dm_f = 0.020;                    % mass flow rate of air (kg/sec) - 0.018
dm_w = 0.017;                    % mass flow rate of water (kg/sec) - 0.015
[t1, x1] = ode23(@(t, x) odeq(t, x, Ta, G, V_w, L, W, dm_f, dm_w), tspan, x0);

T_amb = Ta(floor(t1) + 1);
T_fin = T_amb;
T_win = T_amb;

A_mod = L*W;                          % Module Area
eta_ref = 0.15;                       % Reference Efficiency
beta_p = 0.0045;                      % Temperature coefficient
Tcel_ref = 25;
I_ref = 1000;
delta = 0.052;                        % Solar radiation coefficient
C_f = 1007;                           % Heat capacity (J/kgK)
C_w = 4186;

Two_1 = 2*x1(:,5).' - T_win;
Tfo_1 = 2*x1(:,7).' - T_fin;

Tg_1 = x1(:,1).';
Tc_1 = x1(:,2).';
Tt_1 = x1(:,3).';
Trh_1 = x1(:,4).';
Tw_1 = x1(:,5).';
Trl_1 = x1(:,6).';
Tf_1 = x1(:,7).';
Tp_1 = x1(:,8).';
Ti_1 = x1(:,9).';
groupSize = 3600; % Number of elements in each group
% Calculate the number of groups
numGroups = floor(numel(Two_1)/groupSize);
% Calculate the average of each column (group)
Two = mean(reshape(Two_1(1:numGroups*groupSize), groupSize, numGroups));
Tfo = mean(reshape(Tfo_1(1:numGroups*groupSize), groupSize, numGroups));
Tg = mean(reshape(Tg_1(1:numGroups*groupSize), groupSize, numGroups));
Tc = mean(reshape(Tc_1(1:numGroups*groupSize), groupSize, numGroups));
Tt = mean(reshape(Tt_1(1:numGroups*groupSize), groupSize, numGroups));
Trh = mean(reshape(Trh_1(1:numGroups*groupSize), groupSize, numGroups));
Tw = mean(reshape(Tw_1(1:numGroups*groupSize), groupSize, numGroups));
Trl = mean(reshape(Trl_1(1:numGroups*groupSize), groupSize, numGroups));
Tf = mean(reshape(Tf_1(1:numGroups*groupSize), groupSize, numGroups));
Tp = mean(reshape(Tp_1(1:numGroups*groupSize), groupSize, numGroups));
Ti = mean(reshape(Ti_1(1:numGroups*groupSize), groupSize, numGroups));
T_fin_h = mean(reshape(T_fin(1:numGroups*groupSize), groupSize, numGroups));
T_win_h = mean(reshape(T_win(1:numGroups*groupSize), groupSize, numGroups));

%% Efficiency
% Electrical
eta_E = eta_ref.*(1 - (beta_p.*(Tc - Tcel_ref))+(delta*log(Ia/I_ref)));
eta_E(Ia<=1) = 0;
% Total Thermal
eta_tha = (dm_f.*C_f.*(Tfo - T_fin_h) + dm_w.*C_w.*(Two - T_win_h))./(A_mod.*Ia);
eta_tha(eta_tha<0) = 0;
eta_tha(eta_tha>0.65) = 0.65;
eta_tha(Ia<=1) = 0;
% Water
eta_tha_w = (dm_w.*C_w.*(Two - T_win_h))./(A_mod.*Ia);
eta_tha_w(eta_tha_w<0) = 0;
eta_tha_w(eta_tha_w>0.6) = 0.6;
eta_tha_w(Ia<=1) = 0;
% Air
eta_tha_f = (dm_f.*C_f.*(Tfo - T_fin_h))./(A_mod.*Ia);
eta_tha_f(eta_tha_f<0) = 0;
eta_tha_f(eta_tha_f>0.1) = 0.1;
eta_tha_f(Ia<=1) = 0;

Quele = Ia.*eta_E;
Quth = Ia.*eta_tha;

pvt_collector_output.y(:, 1) = Tg;         % Store Tg values in the first column of y
pvt_collector_output.y(:, 2) = Tc;         % Store Tc values in the second column of y
pvt_collector_output.y(:, 3) = Tt;         % Store Tt values in the third column of y
pvt_collector_output.y(:, 4) = Trh;        % Store Trh values in the fourth column of y
pvt_collector_output.y(:, 5) = Tw;         % Store Tw values in the fifth column of y
pvt_collector_output.y(:, 6) = Trl;        % Store Trl values in the sixth column of y
pvt_collector_output.y(:, 7) = Tf;         % Store Tf values in the seventh column of y
pvt_collector_output.y(:, 8) = Tp;         % Store Tp values in the eight column of y
pvt_collector_output.y(:, 9) = Ti;         % Store Ti values in the ninth column of y
pvt_collector_output.y(:, 10) = Tfo;       % Store Tfo values in the tenth column of y
pvt_collector_output.y(:, 11) = Two;       % Store Two values in the eleventh column of y
pvt_collector_output.y(:, 12) = eta_tha;   % Store eta_tha values in the twelvth column of y
pvt_collector_output.y(:, 13) = eta_E;     % Store eta_E values in the thirteenth column of y

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
plot(xx,T_fin_h,'DisplayName','T_{a}')
hold on
plot(xx, Tg,'DisplayName','T_{gla}')
plot(xx, Tc,'DisplayName','T_{cel}')
plot(xx, Tt,'DisplayName','T_{ted}')
plot(xx, Trh,'DisplayName','T_{rh}')
plot(xx, Tw,'DisplayName','T_{w}')
plot(xx, Trl,'DisplayName','T_{rl}')
plot(xx, Tf,'DisplayName','T_{f}')
plot(xx, Tp,'DisplayName','T_{p}')
plot(xx, Ti,'DisplayName','T_{ins}')
plot(xx, Tfo,'DisplayName','T_{fo}')
plot(xx, Two,'DisplayName','T_{wo}')
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
    title('Power Output Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

    f = figure(1017);
plot(xx,eta_tha_w,'DisplayName','\eta_{w}')
ylabel('Efficiency of Water');
yyaxis right
plot(xx, eta_tha_f,'DisplayName','\eta_{f}')
ylabel('Efficiency of Air');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Efficiency of Water & Air');
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

function X = odeq(t, x, Ta, G, V_w, L, W, dm_f, dm_w)
t2 = floor(t + 1);
Is = G(t2);
T_amb = Ta(t2);
V_w2 = V_w(t2);
T_fin = T_amb;
T_win = T_amb;

w = W*0.85;
d = 0.05;
A_mod = L*W;                     % Module Area
A_c = L*w;                       % Horizontal area of air duct and collector area
A_w = A_mod*0.55;
alpha_gla = 0.1;                 % Absorptivity of glass cover
l_gla = 0.0032;  
rou_gla = 2200;
C_gla = 670;               
M_gla = rou_gla*l_gla*A_mod;
sigma = 5.67e-8;
em_gla = 0.93;                   % Emissivity of glass cover
l_gla = 0.0032;                  % Thickness of glass cover (m)
l_cel = 0.0003;                  % Thickness of solar cell (m)
lamda_gla = 1.1;                 % thermal conductivity of glass cover (W/mK)
lamda_cel = 148;                 % thermal conductivity of solar cell (W/mK)
rou_cel = 2330;
M_cel = rou_cel*l_cel*A_c;
C_cel = 900;
tau_gla = 0.90;                  % Transmisivity of glass cover
alpha_cel = 0.92;                % Absorptivity of solar cell
beta = 0.88;                     % Packing factor
l_ted = 0.00025;                 % Thickness of Tedlar (m)
rou_ted = 1200;         
M_ted = rou_ted*A_mod*l_ted;  
C_ted = 1200;
alpha_ted = 0.6;                 % Absorptivity of tedlar
PC = 0.85;                       % Percent of collector occupied by channel
lamda_ted = 0.2;                 % Thermal conductivity of Tedlar (W/mK)
D_h = 2*(w*d)/(w+d);             % Hydraulic Diameter
mu_f = 1.849e-5;                 % Dynamic viscosity of fluid (kg/m.sec)
lamda_f = 0.026;
lamda_ins = 0.034;               % Thermal conductivity of insulator (W/mK)
l_ins = 0.05;
C_f = 1007;
rou_ins = 20;
V = 1.5;
em_abs = 0.25;
lamda_abs = 160;                 % for Aluminum
l_abs = 0.001;
rho_abs = 2710;
C_abs = 887;
M_abs = rho_abs*A_c*l_abs;
l_w = 0.0015;
rou_w = 998;
C_w = 4186;
M_w = rou_w*l_w*A_w;
M_ins = rou_ins*l_ins*A_mod;
C_ins = 670;
beta_p = 0.0045;                   % Temperature coefficient
delta = 0.052;                     % Solar radiation coefficient
eta_ref = 0.15;                    % Reference Efficiency
Tcel_ref = 25;
I_ref = 1000;
rou = 1.184;
M_f = rou*d*A_c;
% em_ins = 0.6;                    % Emissivity of insulator
% em_ted = 0.85;                   % Emissivity of tedlar
% d_i = 0.008;
% lamda_w = 0.6; 
% D_h_w = d_i;
% mu_w = 0.001308;

Tsk = 0.0552* (T_amb+273.15)^1.5;  % Temperature of Sky
Tsky = Tsk - 273.15;
hva(V_w2<5) = 5.7+3.8*V_w2;        % Convective heat coefficient due to wind
hva(V_w2>=5) = 6.47+V_w2^0.78;

Re = (rou*V*D_h/mu_f);             % Reynolds number
Pr = (mu_f*C_f)/lamda_f;           % Prandtl number
Nu = 0.023*(Re^0.8)*(Pr^0.4);      % Nusselt number

% % % % Dittus–Boelter equation (Forced convection in turbulent pipe flow)
% % % % The Dittus–Boelter equation is valid for:0.6≤Pr≤160, Re≳10000
% % % Re_w = (rou_w*V*D_h_w/mu_w);       % Reynolds number
% % % Pr_w = (mu_w*C_w)/lamda_w;         % Prandtl number
% % % Nu_w = 0.023*(Re_w^0.8)*(Pr_w^0.4);% Nusselt number (0.4: for the fluid being heated)
% % % hvrhmw = Nu_w*((lamda_w)/D_h_w);   % heat transfer coefficient b/w absorber and fluid is fixed equal to 65.1 W/(m2·K)

hcgmc = 1/((l_gla/lamda_gla) + (l_cel/lamda_cel));
hccmg = hcgmc;
hccmt = 1/((l_cel/lamda_cel) + (l_ted/lamda_ted));
hcrhmt = 1/((l_abs/lamda_abs) + (l_ted/lamda_ted));
hcrhmrl = 1/((l_abs/lamda_abs) + (l_abs/lamda_abs));
hvrhmw = 75;
hvrlmw = 10;
hvrlmf = Nu*((lamda_f)/D_h);
hvpmf = hvrlmf;
hcpmi = 1/((l_abs/lamda_abs) + (l_ins/lamda_ins));

%% for glass : x(1)
hrgms = sigma*em_gla*(x(1)^2 + Tsky^2)*(x(1) + Tsky);
x1 = ((alpha_gla*Is) + (hrgms*(Tsky - x(1))) + (hva*(T_amb - x(1))) - ...
    (hcgmc*(x(1) - x(2))))*(A_mod/(M_gla*C_gla));

%% for PV solar cell : x(2)
Quele = Is*A_mod*eta_ref*(1 - (beta_p*(x(2) - Tcel_ref))+(delta*log(Is/I_ref)));
Quele(Is<=0.5) = 0;
x2 = ((tau_gla*alpha_cel*Is*beta) - (hccmg*(x(2) - x(1)) ) - ...
    (hccmt*(x(2) - x(3))))*(A_mod/(M_cel*C_cel)) - (Quele/(M_cel*C_cel));

%% for tedlar : x(3)
x3 = ((tau_gla*alpha_ted*Is*(1 - beta)) + (hccmt*(x(2) - x(3))) - ...
    (hcrhmt*(x(3) - x(4))))*A_mod/(M_ted*C_ted);

%% for higher absorber : x(4)
hrrlmrh = sigma*((x(4) + x(6))*(x(4)^2 + x(6)^2))/((1/em_abs)+(1/em_abs) - 1);
x4 = ((hcrhmt*(x(3) - x(4))) + ((1 - PC)*hcrhmrl*(x(6) - x(4))) + ...
    (PC*hrrlmrh*(x(6) - x(4))))*A_mod/(M_abs*C_abs) + (hvrhmw*(x(5) - x(4)))*A_w/(M_abs*C_abs);

%% for water : x(5)
x5 = ((hvrhmw*(x(4) - x(5))) + (hvrlmw*(x(6) - x(5))))*A_w/(M_w*C_w) - (dm_w*C_w*2*(x(5) - ...
    T_win))*1/(M_w*C_w);

%% for lower absorber : x(6)
hrrlmp = sigma*((x(8) + x(6))*(x(8)^2 + x(6)^2))/((1/em_abs)+(1/em_abs) - 1);
x6 = (((1 - PC)*hcrhmrl*(x(4) - x(6))) + PC*(hrrlmrh*(x(4) - x(6))) + ...
    hrrlmp*(x(8) - x(6)))*A_mod/(M_abs*C_abs) + (hvrlmw*(x(5) - x(6)))*A_w/(M_abs*C_abs) + hvrlmf*(x(7) - ...
    x(6))*A_c/(M_abs*C_abs);

%% for air : x(7)
x7 = ((hvrlmf*(x(6) - x(7))) + (hvpmf*(x(8) - x(7))))*A_c/(M_f*C_f) - (dm_f*C_f*2*(x(7) - ...
    T_fin))*1/(M_f*C_f);

%% for plate+fins : x(8)
x8 = ((hrrlmp*(x(6) - x(8))) + hcpmi*(x(9) - x(8)))*A_mod/(M_abs*C_abs) + (hvpmf*(x(7) - ...
    x(8)))*A_c/(M_abs*C_abs);

%% for insulator : x(9)
x9 = ((hcpmi*(x(8) - x(9))) - (hva*(x(9) - T_amb)))*A_mod/(M_ins*C_ins);
X  = [x1; x2; x3; x4; x5; x6; x7; x8; x9];

end
