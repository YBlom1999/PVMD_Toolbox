function [pvt_collector_output] = bifluid_air_water_tank_redtemp(~, WEATHER_output)
%This function calculates the thermal & PV yield of tank integrated with bi-fluid PVT
% Bi-fluid PV-T with tank (estimated coefficients from real model)
%
% Parameters
% ----------
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% pvt_collector_output : struct
%   Simulation results of the PVT module
% 
% Developed by ZUA.

I = mean(WEATHER_output.Irr,2)';
Tam = [WEATHER_output.ambient_temperature]';
Vw = [WEATHER_output.wind_speed]';

dm_l1 = [0.000833 0.000555 0.0004166 0.000277 0.000222 0.00048 0.0005 0.001527 0.003055 0.0025 0.0020833 0.001944 0.001666 0.0015277 0.00125 0.001 0.0011 0.001388 0.001944 0.0025 0.0023611 0.0020833 0.001666 0.00145833];
dm_l = repmat(dm_l1, [1, 366]);

idx = find(I~=0);                 % find the indices of non-zero elements in I
if length(idx) > length(I)        % if there are non-zero elements in I
    idx = idx(1:length(I));       % select the non-zero elements
end

I1 = I(idx);                      % select the corresponding elements from I

% Second Order Thermal Efficiency Equation (water) unique coefficients
a1a = 8.6667;                     % unique coefficient W/(m2K)                                          
etao_tha_w = 0.5;                 % unique coefficient
C_w = 4186;                       % specific heat of water (J/kg°C) 

% Second Order Thermal Efficiency Equation (air) unique coefficients
a2b = 2.3334;                     % unique coefficient W/(m2K)                                          
etao_tha_f = 0.13;                % unique coefficient
C_f = 1007;                       % specific heat of air (J/kg°C) 

% Electrical Efficiency Equation unique coefficients
e1a = 0.433;                     % unique coefficient W/(m2K) 
etao_ele = 0.15;                 % unique coefficient

DemT= 60;
IntT = Tam(1);

% Calculate surface area of tank
r = 0.4;
h = 0.6;
rho_wat = 1000;                      % water density (kg/m^3)
A_cylinder = 2*pi*r*h;               % Surface area of the cylindrical portion
A_top = pi*r^2;                      % Surface area of the top
A_bottom = pi*r^2;                   % Surface area of the bottom
A_t = A_cylinder+A_top+A_bottom;     % Total surface area
% Calculate volume of tank
vt = pi*r^2*h;                       % Volume of the cylindrical portion
mt= vt*rho_wat;                      % Mass of storage tank (kg) => 1m3=1000kg
L_c = h;                             % characterstic length (m)

Twin = 12*ones(size(Tam));
A_g = 4.33;                      % Gross Collector Area (m2) => 3.28*1.32
dm_w = 0.015;                    % Water Mass Flowrate (kg/s)
dm_f = 0.10;                     % Air Mass Flowrate (kg/s)
T0 = IntT;                       % Initial Temperature of the fluid entering the collector (°C)
Tt0 = IntT;                      % Initial Temperature of the fluid inside tank (°C)
% T_ref = 25;                    % Reference temperature - STC (°C)
% be = 0.0045;                   % Temperature coefficient
d_insulation = 0.08;             % insulation thickness (m)
k_insulation = 0.034;            % insulation thermal conductivity (W/m.K)
rho_air = 1.204;                 % density of air (kg/m^3)
mu_air = 1.82e-5;                % dynamic viscosity of air (Pa.s)
k_air = 0.0257;                  % thermal conductivity of air (W/m.K)
Pr = 0.707;                      % Prandtl number for air
g = 9.81;                        % gravitational acceleration (m/s2)
beta = 257e-6;                   % thermal expansion coefficient (1/T)
d_o = 0.025;                     % outer diameter of tube (m)
v = 8.917e-7;                    % kinematic viscosity (m2/s)
alpha_d = 1.455e-7;              % thermal diffusivity (m2/s)
k = 0.62;                        % thermal conductivity (W/mK)
C = 0.52;                        % constant
n = 0.25;                        % constant (1/4)

% preallocate arrays
TT = zeros(length(Tam), 1);
Tt = zeros(size(Tam));            % storage tank temperature (°C)
T_i = zeros(size(Tam));           % inlet temperature water (°C)
T_fi = zeros(size(Tam));          % inlet temperature air (°C)
eta_tha_w = zeros(size(Tam));
eta_tha_f = zeros(size(Tam));
eta_tha = zeros(size(Tam));
eta_E = zeros(size(Tam));
Tc = zeros(size(Tam));
Qa = zeros(size(Tam));
Two = zeros(size(Tam));
Tfo = zeros(size(Tam));
qctota = zeros(size(Tam));
qloss = zeros(size(Tam));
qload = zeros(size(Tam));
q_hl = zeros(size(Tam));
q_lcmax = zeros(size(Tam));
q_ls = zeros(size(Tam));
Q_aux = zeros(size(Tam));

qaux_w = zeros(size(Tam));
qaux_f = zeros(size(Tam));
Thx_out = zeros(size(Tam));
Tmean = zeros(size(Tam));
tmp = zeros(size(Tam));
tmp_f = zeros(size(Tam));

Diff = zeros(size(Tam));
R_a = zeros(size(Tam));
h_e = zeros(size(Tam));
U_A = zeros(size(Tam));
NTU = zeros(size(Tam));
eta = zeros(size(Tam));
Re = zeros(size(Tam));
Nu = zeros(size(Tam));
h_wind = zeros(size(Tam));
U_t = zeros(size(Tam));

Tt(1) = Tt0;
T_i(1) = T0;

for i = 1:length(Tam)

Re(i) = rho_air*Vw(i)*L_c/mu_air;        % Reynold number
Nu(i) = 0.664*sqrt(Re(i))*Pr^(1/3);      % Nusselt number
h_wind(i) = Nu(i)*k_air/L_c;             % convective heat coefficient due to wind
% Overall heat transfer coefficient of tank
U_t(i) = 1/(1/h_wind(i) + d_insulation/k_insulation);
% U_c = U_t;

% calculate the value that is repeatedly used in the equation for eta_tha
tmp(i) = (T_i(i) - Tam(i)) ./ I(i);
T_fi(i) = Tam(i) + 0.25;
tmp_f(i) = (T_fi(i) - Tam(i)) ./ I(i);
% Second Order Thermal Efficiency Equation
eta_tha_w(i) = etao_tha_w - a1a * tmp(i);%- a2a * tmp(i).^2;
% Second Order Thermal Efficiency Equation (air)
eta_tha_f(i) = etao_tha_f - a2b * tmp_f(i);
% avoid negative values
eta_tha_w(i) = max(eta_tha_w(i), 0);
% avoid values greater than one
eta_tha_w(i) = min(eta_tha_w(i), 1);

% avoid values greater than one
eta_tha_f(i) = min(eta_tha_f(i), 1);
% avoid negative values
eta_tha_f(i) = max(eta_tha_f(i), 0);

% total thermal efficiency
eta_tha(i) = eta_tha_w(i) + eta_tha_f(i);

% Electrical Efficiency Equation
eta_E(i) = etao_ele - e1a * tmp(i);
% avoid negative values
eta_E(i) = max(eta_E(i), 0);
% avoid values greater than one
eta_E(i) = min(eta_E(i), 0.1286);

% fluid outlet temperature
Two(i) = T_i(i) + ((eta_tha_w(i)*A_g*I(i)))/(dm_w*C_w);

% fluid outlet temperature
Tfo(i) = (T_fi(i)) + ((eta_tha_f(i)*A_g*I(i)))/(dm_f*C_f);

% cell temperature
NOCT = 42;                       % typical module at 48°C (best module operated at a NOCT of 33°C, the worst at 58°C)
S = I*0.1;                       % insolation in mW/cm^2
Tc = Tam + S*(NOCT - 20)/80;  
Diff(i) = Two(i) - Tt(i);
if Diff(i)<0
    Diff(i) = 0.05;
end
R_a(i) = (g*beta*(Diff(i))*d_o^3)/(v*alpha_d); % Rayleigh number
h_e(i) = (k*C*R_a(i).^n)/d_o;                       % heat transfer coefficient (W/m2K)
A_e = pi*d_o*L_c;                                   % heat transfer area
% U_A = 1/((1/(h_e*A_e)) + (1/(h_i*A_i)));          % overall heat transfer coefficient
U_A(i) = h_e(i)*A_e;                                % overall heat transfer coefficient
NTU(i) = U_A(i)/(dm_w*C_w);                         % number of transfer units
eta(i) = 1 - exp(-NTU(i));

% Instantaneous total power at normal incidence
Qa(i) = (eta_tha_w(i) + eta_tha_f(i))*A_g*I(i);
% % Instantaneous electrical power at normal incidence
% Qe = eta_E .* A_g .* G;
if i < length(Tam) 
% Collector to tank
qctota(i) = eta(i)*dm_w*C_w*(Two(i) - Tt(i));
% Tank loss
qloss(i) = U_t(i)*A_t*(Tt(i) - Tam(i));
% DHW load
qload(i) = dm_l(i)*C_w*(Tt(i) - Twin(i)); 
UA_l = 36;                                          % space loss coefficient & area product - W/*C (usually 0.2 to 0.5 W/m2*C) 
T_r = 21;                                           % desired air room temperature
eta_l = 0.7;                                        % heat exchange effectiveness
dm_a = 1.1;                                         % flow rate
% Space heating load if considered from tank - (desired room temperature = 21*C)
q_hl(i) = UA_l*(T_r - Tam(i));
q_hl(i) = max(q_hl(i), 0); % Take only positive values by setting negative values to 0
% Max rate of heat transfer
q_lcmax(i) = eta_l*dm_a*C_f*(Tt(i) - T_r);
q_lcmax(i) = max(q_lcmax(i), 0); % Take only positive values by setting negative values to 0
% space load (minimum of the above two - q_hl & q_lcmax)
q_ls(i) = min(q_hl(i), q_lcmax(i));

% tank inlet fluid is already warmer than tank fluid
if Twin(i)>Tt(i)
qload(i) = 0;
end
% qloss
if qloss(i)<0
    qloss(i) = 0;
end
% No heat addition to tank from collector
if Tt(i)>Two(i)
qctota(i) = 0;
end
% inorder to cover hot water demand, auxilary heater is used and addition heat necessary is,
qaux_w(i) = dm_l(i)*C_w*(DemT - Tt(i));
if qaux_w(i)<0
    qaux_w(i) = 0;
end
% if no tank: space heating
qaux_f(i) = dm_f*C_f*(T_r - Tfo(i));
qaux_f(i) = max(qaux_f(i), 0); % Take only positive values by setting negative values to 0

% Total auxilary energy required to cover DHW heating and space loads from tank (positive values only)
Q_aux(i) = q_hl(i) + qaux_w(i) - qloss(i) - q_ls(i);
delta_T = 3600;
Tt(i+1) = Tt(i) + delta_T*(qctota(i) - qloss(i) - qload(i))/(mt*C_w);          % heat loss from tank to environement: (U_t*A_t*th*(Tt(i) - Ta(i)))
end
Thx_out(i) = Two(i) - eta(i)*(Two(i) - Tt(i));
% As it is closed loop
T_i(i+1) = Thx_out(i);

% mean temperature
Tmean(i) = (T_i(i) + Two(i))/2;
% reduced temperature
TT(i) = (Tmean(i) - Tam(i))/I(i);
end

%% Fraction of thermal energy demand covered by PVT system
T_dem = DemT*ones(size(Tt));     
dm_lbb = 0.0014236;
% if Twin<Tt
TA = sum(dm_lbb.*C_w.*(Tt - Twin));       % Load
% else
% TA = 0;
% end
TB = sum(dm_lbb.*C_w.*(T_dem - Twin));    % Demand
F_th = 100.*(TA./TB);
% % % % disp(['Heat remove to cover DHW demand:', num2str(TA/1000), 'kW']);
% % % % disp(['Annual DHW demand:', num2str(TB/1000), 'kW']);
disp(['Fraction of thermal energy demand covered for DHW:', num2str(F_th), '%']);

%%
% Total solar irradiance received by the PVT system in kWh/m^2/year
total_I_sun = sum(mean(WEATHER_output.Irr,2))/1000; 
% Average thermal efficiency
etaSum = sum(eta_tha)/length(I1);
% Average electrical efficiency
etaE = sum(eta_E)/length(I1);

%% change if kWh or kWh/m2
% Yearly total thermal output in kWh/m2
total_thermal_output = etaSum*total_I_sun*A_g/A_g;
% Yearly total electrical output in kWh/m2
total_electrical_output = etaE*total_I_sun*A_g/A_g;

pvt_collector_output.eta_tha=eta_tha;
pvt_collector_output.etaSum=etaSum;
pvt_collector_output.eta_E=eta_E;
pvt_collector_output.etaE=etaE;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.Tt=Tt;
pvt_collector_output.Tmean=Tmean;
pvt_collector_output.total_thermal_output=total_thermal_output;
pvt_collector_output.total_electrical_output=total_electrical_output;

%% Plot results
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
plot(xx,eta_tha_w,'DisplayName','\eta_{th,w}')
hold on 
plot(xx,eta_tha_f,'DisplayName','\eta_{th,f}')
plot(xx,eta_tha,'DisplayName','\eta_{tot}')

ylabel('temperature [^oC]')
xlabel('month')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Temperature Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

f = figure(100);
plot(xx,Thx_out,'DisplayName','T_{hx,o}')
hold on
plot(xx,Tmean,'DisplayName','T_{m}')
plot(xx,Two,'DisplayName','T_{w,o}')
plot(xx,Tfo,'DisplayName','T_{f,o}')
plot(xx,Tt,'DisplayName','T_{t}')
plot(xx,Tam,'DisplayName','T_{am}')
% plot(xx,Tc,'DisplayName','T_{c}')
ylabel('temperature [^oC]')
xlabel('month')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Temperature Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

f = figure(1021);
plot(xx,qload,'DisplayName','Q_{load}')
ylabel('power output [W]')
yyaxis right
plot(xx,qaux_w, 'DisplayName','Q_{aux}')
ylabel('auxilary heater power [W]')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Load & Auxilary Power Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

    f = figure(1023);
plot(xx,q_hl,'DisplayName','Q_{sp,load}')
ylabel('power output [W]')
yyaxis right
plot(xx,qaux_f, 'DisplayName','Q_{sp,aux}')
ylabel('auxilary heater power [W]')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Load & Auxilary Power Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

    % ELectical power
%     P_el = zeros(size(Tam));
    P_el = eta_E.*I.*A_g;
    f = figure(1024);
plot(xx,P_el,'DisplayName','P_{el}')

ylabel('Electrical power [W]')
xlabel('month')
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Electrical power');
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
