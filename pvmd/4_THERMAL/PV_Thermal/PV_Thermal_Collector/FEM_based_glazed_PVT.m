function [pvt_collector_output] = FEM_based_glazed_PVT(~,MODULE_output, WEATHER_output)
% Calculates the thermal and PV yield of PVT collector using finite element method (FEM)
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
% Developed by J & ZUA.
%
Ia = mean(WEATHER_output.Irr,2)';
Tam = [WEATHER_output.ambient_temperature]';
L = MODULE_output.ML;
W = MODULE_output.MW;  
dm_w = 0.02;                          % mass flow rate (kg/sec)
Ia(Ia<=0.05) = 0.05;                  % so that efficiency doesn't give NaN values

idx = find(Ia>0.05);                  % find the indices of non-zero elements in Ia
if length(idx) > length(Ia)           % if there are non-zero elements in Ia
    idx = idx(1:length(Ia));          % select the non-zero elements
end
I1 = Ia(idx);                         % length when there is solar radiation > 0.05  

%% variables
Start_temp = 20.15;
% Output vectors
Two = zeros(1,length(Ia)); 
Tc = zeros(1,length(Ia)); 
eta_E = zeros(1,length(Ia));
eta_tha = zeros(1,length(Ia));
eta_tot = zeros(1,length(Ia));
Quele = zeros(1,length(Ia));
Quth = zeros(1,length(Ia));
Power_loss_day = zeros(1,length(Ia));
    
% Calculations per timestep
for hour = 1:1:length(Ia)

Irr = Ia(hour);                         % [W/m2]
Twin = 20 + 273;                        % [°K]
 
    if Irr < 100
      Two(hour) = Twin - 273;           % [°C] 
      Tc(hour) = Twin - 273;            % [°C]
      eta_E(hour) = 0;
      eta_tha(hour) = 0;
      eta_tot(hour) = 0;
      Quele(hour) = 0;                  % [W]
      Quth(hour) = 0;                   % [W]
      Power_loss_day(hour) = 0;         % [W]
      Tc_1 = Start_temp;                % [°C] 
    else
      th_pipe = 0.004;                  % pipe thickness [m] 
      space_pipe = 0.095;               % pipe spacing [m] 
      N_pipe = 15;                      % number of pipes
      th_insX = 0.02;                   % insulation wall [m]
      th_insY = 0.02;                   % insulation wall [m]
      Phi = 31;                         % module inclination [°]
      Tc_1 = Start_temp;                % [°C] 
      T_PV = Tc_1 + 273;                % [°K]       
      T_ambient = Tam(hour) + 273;      % [°K]
    
      Xn_max = 800;                     % number of points in X direction
      Yn_max = 25;                      % number of points in y direction
      Zn_max = 10;                      % number of points in Z direction in middle section
      totalsteps = 300;                 % number of calculation cycles to guess steady state for middle section
      stepc = 600;                      % circle step size, higher numbers lead to smoother curves

% th_glassc = 0.0032;                   % thickness glass cover [m]
% th_airgap = 0.02;                     % thickness air gap cover [m]
% th_glassPV = 0.0015;                  % thickness glass cover [m]
% th_EVAt = 0.0015;                     % thickness PV layer [m]
% th_PV = 0.00035;                      % thickness PV layer [m]
th_EVAb = 0.0015;                       % thickness PV layer [m]
th_TED = 0.001;                         % thickness PV layer [m]
th_ad = 0.001;                          % thickness PV layer [m]
th_cu = 0.002;                          % thickness copper layer [m]

plotY = th_TED + th_ad + th_cu+th_insY; % height of collector, Y direction [m]       
rad = 0.005;                            % inner radius of the pipe [m]
eta_panel = 0.72;
sigma = 5.67*10^-8;
I_ref = 1000; 
Tcel_ref = 25;
N_cells = 68;
k_nu = -0.0029;
Voc = 48.3080;
Isc = 6.39;
FF = 0.79;
kb = 8.61*10^-5; 
% A_mod = L*W;
% A_pipe = pi*(rad^2);                 % inner surface of the pipe [m^2] 
% eta_sky = 0.68;         
% Vw = 1; 

% Definitions
dX = L/(Xn_max);                       % stepsize X-direction
dY = plotY/(Yn_max);                   % stepsize Y-direction
dZ = W/(Zn_max);                       % stepsize Z-direction
% Vm = dX*dY*dZ;                       % volume of block middle

Xn_ins = round(th_insX/dX);            % number of points of insulation material in X direction
Yn_ins = round(th_insY/dY);            % number of points of insulation material in Y direction
Yn_TED = round(th_TED/dY);             % thickness PV layer [m]
Yn_ad = round(th_ad/dY);               % thickness PV layer [m]
Yn_cu = round(th_cu/dY);               % number of points of insulation material in Y direction
Xn_tot = Xn_max+2*Xn_ins+2;            % total number of points in X direction
Yn_tot = Yn_max+2;                     % total number of points in Y direction
Zn_tot = Zn_max+2;                     % total number of points in Z direction
% Yn_EVAb = round(th_EVAb/dY);         % thickness PV layer [m]

% Lower boundry test
th_EVAbtestY = th_EVAb/dY;             % thickness PV layer [m]
th_TEDtestY = th_TED/dY;               % thickness PV layer [m]
th_adtestY = th_ad/dY;                 % thickness PV layer [m]
radtestX = rad/dX;
radtestY = rad/dY;
th_pipetestX = th_pipe/dX;
th_pipetestY = th_pipe/dY;
th_instestX = th_insX/dX;
th_instestY = th_insY/dY;

if radtestX < 1
    disp(radtestX)
    error('Pipe radius less then 1 point in X direction')
elseif th_pipetestX < 1
    disp(th_pipetestX)
    error('Pipe thickness less then 1 point in X direction')
elseif th_instestX < 1
    disp(th_instestX)
    error('Insulation thickness less then 1 point in X direction')
elseif radtestY < 1
    disp(radtestY)
    error('Pipe radius then 1 point in Y direction')
elseif th_pipetestY < 1
    disp(th_pipetestY)
    error('Pipe thickness less then 1 point in Y direction')
elseif th_instestY < 1
    disp(th_instestY)
    error('Insulation thickness less then 1 point in Y direction')
elseif th_EVAbtestY < 1
    disp(th_EVAbtestY)
    error('bottom EVA thickness less then 1 point in Y direction')
elseif th_TEDtestY < 1
    disp(th_TEDtestY)
    error('TED thickness less then 1 point in Y direction')
elseif th_adtestY < 1
    disp(th_adtestY)
    error('adhesive thickness less then 1 point in Y direction')
end

%% Material properties
K_cu = 389;                            % [W/mK] thermal conductivity copper
K_al = 237;                            % [W/mK] thermal conductivity aluminum
K_PV = 148;                            % [W/mK]  
K_EVAt = 0.23;                         % [W/mK]  
K_EVAb = 0.23;                         % [W/mK]  
K_TED = 0.07;                          % [W/mK]  
K_ad = 1.5;                            % [W/mK]  
K_glassc = 0.6;                        % [W/mK] 
eta_glassc = 0.08;                     % transmittence
K_glassPV = 1.1;                       % [W/mK] 
K_ins = 0.034;                         % [W/mK]
K_fluid = 0.6;                         % [W/mK] 
C_w = 4182;                            % [J/KgK] 
K_air = 0.0261;                        % [W/mK] 
hc = 0.5*(48/11)*K_fluid/rad;          % convection coefficient water

% rho_cu = 8900;                       % [kg/m^3] 
% Cp_cu = 386;                         % [J/KgK]
% rho_al = 2702;                       % [kg/m^3] 
% Cp_al = 880;                         % [J/KgK]
% rho_PV = 2330;                       % [kg/m^3]  
% Cp_PV = 760;                         % [J/KgK] 
% rho_EVAt = 921;                      % [kg/m^3]  
% Cp_EVAt = 2300;                      % [J/KgK] 
% rho_EVAb = 921;                      % [kg/m^3]  
% Cp_EVAb = 2300;                      % [J/KgK] 
% rho_TED = 1200;                      % [kg/m^3]  
% Cp_TED = 1250;                       % [J/KgK] 
% rho_ad = 1200;                       % [kg/m^3]  
% Cp_ad = 1250;                        % [J/KgK] 
% rho_glassc = 2200;                   % [kg/m^3]  
% Cp_glassc = 670;                     % [J/KgK] 
% rho_glassPV = 2200;                  % [kg/m^3]  
% Cp_glassPV = 670;                    % [J/KgK] 
% rho_ins = 20;                        % [kg/m^3]  
% Cp_ins = 670;                        % [J/KgK] 
% rho_fluid = 998;                     % [kg/m^3]  
% nu_fluid = 0.001002;                 % [kg/ms]
% rho_air = 1.220;                     % [kg/m^3]  
% Cp_air = 1007;                       % [J/KgK] 
% K_wa = 0.598;                        % [W/mK] thermal conductivity water
% K_si = 380;                          % [W/mK] thermal pure silicon

%% Empty vectors and initial values
% convection
Beta_air = 0.0034;                     % [1/K]
V_air = 15.66*10^-6;                   % [m^2/s]
g = 9.81;                              % [m/s^2]
Pr = 0.69;
% Alpha_air = V_air/Pr;                % [m^2/s]
% Nu_air = 18.43*10^-6;                % [Kg/m*s]
% T1 = zeros(Yn_tot,Xn_tot,Zn_tot);
        
T_boundaryX = zeros(Yn_tot,Xn_tot,Zn_tot);
T_boundaryY = zeros(Yn_tot,Xn_tot,Zn_tot);
T_boundaryZ = zeros(Yn_tot,Xn_tot,Zn_tot);

Q1 = zeros(Yn_tot,Xn_tot,Zn_tot);
Q_air = zeros(Yn_tot,Xn_tot,Zn_tot);
Q_bottom = zeros(1,Xn_tot,Zn_tot);
Materials = strings(Yn_tot,Xn_tot);
Transfer = strings(Yn_tot,Xn_tot,Zn_tot);

T_sky = T_ambient - 20;
T_air = T_ambient; 
T_glassc = T_PV - 0.01;
T_TED = Twin; 
T_ad = Twin;                            
T_cu = Twin;
T_ins = Twin;
% T_EVAb = Twin;
% T_fluid = Twin;
% T_al = Twin;

%% 2D design
while 1
% Middle section
T = [[repmat(T_air,1,1),repmat(T_ins,1,Xn_ins),repmat(repmat(T_PV,1,Xn_max),1,1),repmat(T_ins,1,Xn_ins),repmat(T_air,1,1)];[repmat(T_air,Yn_tot-2,1),repmat(T_ins,Yn_tot-2,Xn_ins),[repmat(repmat(T_TED,1,Xn_max),Yn_TED,1);repmat(repmat(T_ad,1,Xn_max),Yn_ad,1);repmat(repmat(T_cu,1,Xn_max),Yn_cu,1);repmat(repmat(T_ins,1,Xn_max),Yn_ins,1)],repmat(T_ins,Yn_tot-2,Xn_ins),repmat(T_air,Yn_tot-2,1)];repmat(T_air,1,Xn_tot)];
k = [[repmat(K_air,1,1),repmat(K_ins,1,Xn_ins),repmat(repmat(K_PV,1,Xn_max),1,1),repmat(K_ins,1,Xn_ins),repmat(K_air,1,1)];[repmat(K_air,Yn_tot-2,1),repmat(K_ins,Yn_tot-2,Xn_ins),[repmat(repmat(K_TED,1,Xn_max),Yn_TED,1);repmat(repmat(K_ad,1,Xn_max),Yn_ad,1);repmat(repmat(K_cu,1,Xn_max),Yn_cu,1);repmat(repmat(K_ins,1,Xn_max),Yn_ins,1)],repmat(K_ins,Yn_tot-2,Xn_ins),repmat(K_air,Yn_tot-2,1)];repmat(K_air,1,Xn_tot)];

% Pipes
thetas = linspace(0, (2*pi), stepc);
for npipes = 1:1:N_pipe
    dY_rad = round(rad/dY)+round(th_pipe/dY);     
    dX_rad = round(rad/dX)+round(th_pipe/dX);  
    Y_centre = ceil((th_TED + th_ad + th_cu + th_pipe + rad)/dY);
    X_centre = round((Xn_tot)/2 - 0.5*space_pipe*((1-(-1)^(N_pipe+1))/2)/dX + (space_pipe/dX)*(npipes-ceil((N_pipe/2))));    
    T(Y_centre-dY_rad:Y_centre,X_centre-dX_rad:X_centre+dX_rad) = T_cu;
    k(Y_centre-dY_rad:Y_centre,X_centre-dX_rad:X_centre+dX_rad) = K_al;    

    for rx = 0:1:round(rad/dX)
        for ry = 0:1:round(rad/dY)
            circleXr = round(rx * cos(thetas)+ (Xn_tot)/2 - 0.5*space_pipe*((1-(-1)^(N_pipe+1))/2)/dX + (space_pipe/dX)*(npipes-ceil((N_pipe/2))));
            circleYr = round(ry * sin(thetas)+ ceil((th_TED + th_ad + th_cu + th_pipe + rad)/dY));

            for i = 1:1:stepc     
                T(circleYr(i),circleXr(i)) = Twin;
                k(circleYr(i),circleXr(i)) = K_fluid;                
            end
        end
    end
end

for npipesp = 1:1:N_pipe
    for rpx = round(rad/dX)+1:1:round(rad/dX)+round(th_pipe/dX)
        for rpy = round(rad/dY)+1:1:round(rad/dY)+round(th_pipe/dY)
            circleXrp = round(rpx * cos(thetas)+ (Xn_tot)/2 - 0.5*space_pipe*((1-(-1)^(N_pipe+1))/2)/dX + (space_pipe/dX)*(npipesp-ceil((N_pipe/2))));
            circleYrp = round(rpy * sin(thetas)+ ceil((th_TED + th_ad + th_cu + th_pipe + rad)/dY));
            for ip = 1:1:stepc     

                T(circleYrp(ip),circleXrp(ip)) = T_cu;
                k(circleYrp(ip),circleXrp(ip)) = K_al;             
            end 
        end
    end
end

% Material definitions
for mx = 1:Xn_tot
    for my = 1:Yn_tot
        if k(my,mx) == K_cu
            Materials(my,mx) = 'Cu';
        elseif k(my,mx) == K_al
            Materials(my,mx) = 'Al';
        elseif k(my,mx) == K_glassc
            Materials(my,mx) = 'cover';
        elseif k(my,mx) == K_glassPV
            Materials(my,mx) = 'glassPV';  
        elseif k(my,mx) == K_EVAt
            Materials(my,mx) = 'EVAt';    
        elseif k(my,mx) == K_PV
            Materials(my,mx) = 'PV';
        elseif k(my,mx) == K_EVAb
            Materials(my,mx) = 'EVAb';
        elseif k(my,mx) == K_TED
            Materials(my,mx) = 'TED'; 
        elseif k(my,mx) == K_ad
            Materials(my,mx) = 'adhesive';    
        elseif k(my,mx) == K_ins
            Materials(my,mx) = 'Ins';
        elseif k(my,mx) == K_fluid
            Materials(my,mx) = 'Fluid';
        elseif k(my,mx) == K_air
            Materials(my,mx) = 'Air';   
        end   
    end
end

%% 3D arrays
% Module
T3d = cat(3,repmat(T_air,Yn_tot,Xn_tot),repmat(T,[1 1 Zn_max]),repmat(T_air,Yn_tot,Xn_tot));
k3d = cat(3,repmat(K_air(1),Yn_tot,Xn_tot),repmat(k,[1 1 Zn_max]),repmat(K_air(1),Yn_tot,Xn_tot));

%% formulas
while 1           
    Q_air_rad =(eta_panel*eta_glassc/(eta_panel+eta_glassc-eta_panel*eta_glassc))*sigma*(T_PV^4-T_glassc^4); 

    Ra = Pr*((g*Beta_air*(T_PV - T_glassc).*(L^3))/(V_air^2)); 
    if Ra <0
        error('T_PV too low for T_ambient')
    end
    Ra_crit = 10^8;
    Nu = 1+ 1.44*(1 - (1708*((sind(1.8*Phi))^1.6))/(Ra*cosd(Phi)))*(max(0,(1 - 1708/(Ra*cosd(Phi))))) + max(0,(((Ra*cosd(Phi))/(5830))^0.333 - 1)); 
    hc_air = Nu*K_glassPV/L; 
    Q_air_conv = hc_air*L*W.*(T_PV - T_glassc);  

    Loss_radiation = L*W*eta_panel*sigma*(T_glassc^4 - T_sky^4);

    Ra2 = Pr*((g*Beta_air*(T_glassc - T_ambient).*(L^3))/(V_air^2)); 
    if Ra2 <0
        Ra2 = 0;                     % error('T_glassc too low for T_ambient')
    end
    Nu2 = 0.56*((Ra_crit*cosd(Phi))^0.25) + 0.13*((Ra2^0.333) - (Ra_crit^0.333)); 
    hc_sky = Nu2*K_glassc/L;   
    Loss_convection = hc_sky*L*W.*(T_glassc - T_air); 

    Q_Loss = Q_air_rad + Q_air_conv - Loss_radiation - Loss_convection;

    if Q_Loss < 2 && Q_Loss > -2
        break
    end
    if Q_Loss < 0
        T_glassc = T_glassc - 0.01;
    else
        T_glassc = T_glassc + 0.01;
    end
end

    for n = 1:totalsteps 
    T1 = T3d;
    for X = 2:Xn_tot-1
        for Y = 2:Yn_tot-1
            for Z = 2:Zn_tot-1
                % 3D conduction solids
                if (k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_air && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_air && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_air && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_air && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z+1) ~= K_air && k3d(Y,X,Z-1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_air)
                % Isotropic
                if (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T1(Y,X,Z) = (1/(2/(dX^2)+2/(dY^2)+2/(dZ^2)))*((T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond';
                % Single border boundary conduction
                % X =>
                elseif (k3d(Y,X+1,Z) ~= k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X+1,Z)*T3d(Y,X+1,Z))/(k3d(Y,X,Z)+k3d(Y,X+1,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+2/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X-1,Z))/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X =>';
                % X <=
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) ~= k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X-1,Z)*T3d(Y,X-1,Z))/(k3d(Y,X,Z)+k3d(Y,X-1,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+2/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X+1,Z))/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X <=';
                % Y ^
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) ~= k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y-1,X,Z)*T3d(Y-1,X,Z))/(k3d(Y,X,Z)+k3d(Y-1,X,Z));
                    T1(Y,X,Z) = (1/(2/(dX^2)+3/(dY^2)+2/(dZ^2)))*((T3d(Y,X-1,Z)+T3d(Y,X+1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y+1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond Y ^';
                % Y v
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) ~= k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T1(Y,X,Z) = (1/(2/(dX^2)+3/(dY^2)+2/(dZ^2)))*((T3d(Y,X-1,Z)+T3d(Y,X+1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond Y v';
                % Z +
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) ~= k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryZ(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X,Z+1)*T3d(Y,X,Z+1))/(k3d(Y,X,Z)+k3d(Y,X,Z+1));
                    T1(Y,X,Z) = (1/(2/(dX^2)+2/(dY^2)+3/(dZ^2)))*((T3d(Y,X-1,Z)+T3d(Y,X+1,Z))/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(2*T_boundaryZ(Y,X,Z)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond Z +'; 
                % Z -
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) ~= k3d(Y,X,Z))
                    T_boundaryZ(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X,Z-1)*T3d(Y,X,Z-1))/(k3d(Y,X,Z)+k3d(Y,X,Z-1));
                    T1(Y,X,Z) = (1/(2/(dX^2)+2/(dY^2)+3/(dZ^2)))*((T3d(Y,X-1,Z)+T3d(Y,X+1,Z))/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(2*T_boundaryZ(Y,X,Z)+T3d(Y,X,Z+1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond Z -'; 
                % double border boundary conduction 
                % X =>, Y v
                elseif (k3d(Y,X+1,Z) ~= k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) ~= k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X+1,Z)*T3d(Y,X+1,Z))/(k3d(Y,X,Z)+k3d(Y,X+1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+3/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X-1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X =>, Y v';
                % X <=, Y v
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) ~= k3d(Y,X,Z) && k3d(Y+1,X,Z) ~= k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X-1,Z)*T3d(Y,X-1,Z))/(k3d(Y,X,Z)+k3d(Y,X-1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+3/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X+1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X <=, Y v';
                % X =>, Y ^
                elseif (k3d(Y,X+1,Z) ~= k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) ~= k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X+1,Z)*T3d(Y,X+1,Z))/(k3d(Y,X,Z)+k3d(Y,X+1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y-1,X,Z)*T3d(Y-1,X,Z))/(k3d(Y,X,Z)+k3d(Y-1,X,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+3/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X-1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y+1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X =>, Y ^';
                % X <=, Y ^
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) ~= k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) ~= k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X-1,Z)*T3d(Y,X-1,Z))/(k3d(Y,X,Z)+k3d(Y,X-1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y-1,X,Z)*T3d(Y-1,X,Z))/(k3d(Y,X,Z)+k3d(Y-1,X,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+3/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X+1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y+1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X <=, Y ^';
                % Y v, Y ^
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) ~= k3d(Y,X,Z) && k3d(Y-1,X,Z) ~= k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T_boundaryZ(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y-1,X,Z)*T3d(Y-1,X,Z))/(k3d(Y,X,Z)+k3d(Y-1,X,Z));
                    T1(Y,X,Z) = (1/(2/(dX^2)+4/(dY^2)+2/(dZ^2)))*((T3d(Y,X-1,Z)+T3d(Y,X+1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+(2*T_boundaryZ(Y,X,Z)))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond Y v, Y ^';
                % triple border boundary conduction 
                % X =>, Y v, Y ^
                elseif (k3d(Y,X+1,Z) ~= k3d(Y,X,Z) && k3d(Y,X-1,Z) == k3d(Y,X,Z) && k3d(Y+1,X,Z) ~= k3d(Y,X,Z) && k3d(Y-1,X,Z) ~= k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X+1,Z)*T3d(Y,X+1,Z))/(k3d(Y,X,Z)+k3d(Y,X+1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T_boundaryZ(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y-1,X,Z)*T3d(Y-1,X,Z))/(k3d(Y,X,Z)+k3d(Y-1,X,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+4/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X-1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+(2*T_boundaryZ(Y,X,Z)))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X =>, Y v, Y ^';
                % X <=, Y v, Y ^
                elseif (k3d(Y,X+1,Z) == k3d(Y,X,Z) && k3d(Y,X-1,Z) ~= k3d(Y,X,Z) && k3d(Y+1,X,Z) ~= k3d(Y,X,Z) && k3d(Y-1,X,Z) ~= k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X-1,Z)*T3d(Y,X-1,Z))/(k3d(Y,X,Z)+k3d(Y,X-1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T_boundaryZ(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y-1,X,Z)*T3d(Y-1,X,Z))/(k3d(Y,X,Z)+k3d(Y-1,X,Z));
                    T1(Y,X,Z) = (1/(3/(dX^2)+4/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+T3d(Y,X+1,Z))/(dX^2)+(2*T_boundaryY(Y,X,Z)+(2*T_boundaryZ(Y,X,Z)))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X <=, Y v, Y ^';
                % X =>, X <=, Y ^
                elseif (k3d(Y,X+1,Z) ~= k3d(Y,X,Z) && k3d(Y,X-1,Z) ~= k3d(Y,X,Z) && k3d(Y+1,X,Z) == k3d(Y,X,Z) && k3d(Y-1,X,Z) ~= k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X-1,Z)*T3d(Y,X-1,Z))/(k3d(Y,X,Z)+k3d(Y,X-1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T_boundaryZ(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X+1,Z)*T3d(Y,X+1,Z))/(k3d(Y,X,Z)+k3d(Y,X+1,Z));
                    T1(Y,X,Z) = (1/(4/(dX^2)+3/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+(2*T_boundaryZ(Y,X,Z)))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y+1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X =>, X <=, Y ^';
                % X =>, X <=, Y v
                elseif (k3d(Y,X+1,Z) ~= k3d(Y,X,Z) && k3d(Y,X-1,Z) ~= k3d(Y,X,Z) && k3d(Y+1,X,Z) ~= k3d(Y,X,Z) && k3d(Y-1,X,Z) == k3d(Y,X,Z) && k3d(Y,X,Z+1) == k3d(Y,X,Z) && k3d(Y,X,Z-1) == k3d(Y,X,Z))
                    T_boundaryX(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X-1,Z)*T3d(Y,X-1,Z))/(k3d(Y,X,Z)+k3d(Y,X-1,Z));
                    T_boundaryY(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y+1,X,Z)*T3d(Y+1,X,Z))/(k3d(Y,X,Z)+k3d(Y+1,X,Z));
                    T_boundaryZ(Y,X,Z) = (k3d(Y,X,Z)*T3d(Y,X,Z) + k3d(Y,X+1,Z)*T3d(Y,X+1,Z))/(k3d(Y,X,Z)+k3d(Y,X+1,Z));
                    T1(Y,X,Z) = (1/(4/(dX^2)+3/(dY^2)+2/(dZ^2)))*((2*T_boundaryX(Y,X,Z)+(2*T_boundaryZ(Y,X,Z)))/(dX^2)+(2*T_boundaryY(Y,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = 'cond X =>, X <=, Y v';
                else
                    T1(Y,X,Z) = (1/(2/(dX^2)+2/(dY^2)+2/(dZ^2)))*((T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(dZ^2)); 
                    Transfer(Y,X,Z) = '???';
                end
                % Single border convection pipe-fluid
                % X =>
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc/dX))*(hc*T3d(Y,X+1,Z)/dX+(T3d(Y,X-1,Z)/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'X =>';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z));
                % X <=
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc/dX))*(hc*T3d(Y,X-1,Z)/dX+(T3d(Y,X+1,Z)/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'X <=';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z));
                % Y ^
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc/dY))*(hc*T3d(Y-1,X,Z)/dY+(T3d(Y+1,X,Z)/(dY^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'Y ^';
                    Q1(Y,X,Z) = -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z));
                % Y v
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc/dY))*(hc*T3d(Y+1,X,Z)/dY+(T3d(Y-1,X,Z)/(dY^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'Y v';
                    Q1(Y,X,Z) = -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z));
                % Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc/dZ))*(hc*T3d(Y,X,Z+1)/dZ+(T3d(Y,X,Z-1)/(dZ^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'Z +';
                    Q1(Y,X,Z) = -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc/dZ))*(hc*T3d(Y,X,Z-1)/dZ+(T3d(Y,X,Z+1)/(dZ^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'Z -';
                    Q1(Y,X,Z) = -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));
                % Single border air
                % X =>
                elseif k3d(Y,X+1,Z) == K_air
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc_air/dX))*(hc_air*T3d(Y,X+1,Z)/dX+(T3d(Y,X-1,Z)/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'air X =>';
                    Q_air(Y,X,Z) = -hc_air*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z));
                % X <=
                elseif k3d(Y,X-1,Z) == K_air
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc_air/dX))*(hc_air*T3d(Y,X-1,Z)/dX+(T3d(Y,X+1,Z)/(dX^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'air X <=';
                    Q_air(Y,X,Z) = -hc_air*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z));
                % Y ^
                elseif k3d(Y-1,X,Z) == K_air
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc_air/dY))*(hc_air*T3d(Y-1,X,Z)/dY+(T3d(Y+1,X,Z)/(dY^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'air Y ^';
                    Q_air(Y,X,Z) = -hc_air*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z));
                %Y v
                elseif k3d(Y+1,X,Z) == K_air
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc_air/dY))*(hc_air*T3d(Y+1,X,Z)/dY+(T3d(Y-1,X,Z)/(dY^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'air Y v';    
                    Q_air(Y,X,Z) = -hc_air*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z));
                % double border convection pipe-fluid
                % X =>, Y v
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y+1,X,Z)/dY)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y-1,X,Z)/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'X =>, Y v';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z));
                % X <=, Y v
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y+1,X,Z)/dY)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y-1,X,Z)/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'X <=, Y v';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z));
                % X <=, Y ^
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y-1,X,Z)/dY)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y+1,X,Z)/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'X <=, Y ^';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z));
                % X =>, Y ^
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y-1,X,Z)/dY)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y+1,X,Z)/(dY^2)+(T3d(Y,X,Z+1)+T3d(Y,X,Z-1))/(2*dZ^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'X =>, Y ^';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z));
                % Y ^, Z + 
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dY+1/dZ)))*(hc*(T3d(Y-1,X,Z)/dY+T3d(Y,X,Z+1)/dZ)+(T3d(Y+1,X,Z)/(dY^2)+T3d(Y,X,Z-1)/(dZ^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'Y ^, Z +'; 
                    Q1(Y,X,Z) = -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z)) - hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % Y ^, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dY+1/dZ)))*(hc*(T3d(Y-1,X,Z)/dY+T3d(Y,X,Z-1)/dZ)+(T3d(Y+1,X,Z)/(dY^2)+T3d(Y,X,Z+1)/(dZ^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'Y ^, Z -'; 
                    Q1(Y,X,Z) = -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));
                % Y v, Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dY+1/dZ)))*(hc*(T3d(Y+1,X,Z)/dY+T3d(Y,X,Z+1)/dZ)+(T3d(Y-1,X,Z)/(dY^2)+T3d(Y,X,Z-1)/(dZ^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2))*k3d(Y,X,Z)); 
                    Transfer(Y,X,Z) = 'Y v, Z +'; 
                    Q1(Y,X,Z) = -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % Y v, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dY+1/dZ)))*(hc*(T3d(Y+1,X,Z)/dY+T3d(Y,X,Z-1)/dZ)+(T3d(Y-1,X,Z)/(dY^2)+T3d(Y,X,Z+1)/(dZ^2)+(T3d(Y,X+1,Z)+T3d(Y,X-1,Z))/(2*dX^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'Y v, Z -';
                    Q1(Y,X,Z) = -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));
                % X =>, Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dZ)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y,X,Z+1)/dZ)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y,X,Z-1)/(dZ^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2))*k3d(Y,X,Z));   
                    Transfer(Y,X,Z) = 'X =>, Z +'; 
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % X =>, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dZ)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y,X,Z-1)/dZ)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y,X,Z+1)/(dZ^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2))*k3d(Y,X,Z));   
                    Transfer(Y,X,Z) = 'X =>, Z -'; 
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));
                % X <=, Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dZ)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y,X,Z+1)/dZ)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y,X,Z-1)/(dZ^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2))*k3d(Y,X,Z));   
                    Transfer(Y,X,Z) = 'X <=, Z +'; 
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % X <=, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dZ)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y,X,Z-1)/dZ)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y,X,Z+1)/(dZ^2)+(T3d(Y+1,X,Z)+T3d(Y-1,X,Z))/(2*dY^2))*k3d(Y,X,Z));   
                    Transfer(Y,X,Z) = 'X <=, Z -'; 
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));

                % triple border convection pipe-fluid
                % X =>, Y ^, Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y-1,X,Z)/dY+T3d(Y,X,Z+1)/dZ)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y+1,X,Z)/(dY^2)+T3d(Y,X,Z-1)/(dZ^2))*k3d(Y,X,Z));  
                    Transfer(Y,X,Z) = 'X =>, Y ^, Z +';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % X =>, Y ^, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y-1,X,Z)/dY+T3d(Y,X,Z-1)/dZ)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y+1,X,Z)/(dY^2)+T3d(Y,X,Z+1)/(dZ^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'X =>, Y ^, Z -';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));
                % X =>, Y v, Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y+1,X,Z)/dY+T3d(Y,X,Z+1)/dZ)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y-1,X,Z)/(dY^2)+T3d(Y,X,Z-1)/(dZ^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'X =>, Y v, Z +';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % X =>, Y v, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) == K_fluid && k3d(Y,X-1,Z) ~= K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X+1,Z)/dX+T3d(Y+1,X,Z)/dY+T3d(Y,X,Z-1)/dZ)+(T3d(Y,X-1,Z)/(dX^2)+T3d(Y-1,X,Z)/(dY^2)+T3d(Y,X,Z+1)/(dZ^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'X =>, Y v, Z -';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X+1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));
                % X <=, Y ^, Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y-1,X,Z)/dY+T3d(Y,X,Z+1)/dZ)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y+1,X,Z)/(dY^2)+T3d(Y,X,Z-1)/(dZ^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'X <=, Y ^, Z +';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % X <=, Y ^, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) ~= K_fluid && k3d(Y-1,X,Z) == K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y-1,X,Z)/dY+T3d(Y,X,Z-1)/dZ)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y+1,X,Z)/(dY^2)+T3d(Y,X,Z+1)/(dZ^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'X <=, Y ^, Z -';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y-1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));
                % X <=, Y v, Z +
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) == K_fluid && k3d(Y,X,Z-1) ~= K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y+1,X,Z)/dY+T3d(Y,X,Z+1)/dZ)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y-1,X,Z)/(dY^2)+T3d(Y,X,Z-1)/(dZ^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'X <=, Y v, Z +';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z+1) - T3d(Y,X,Z));
                % X <=, Y v, Z -
                elseif (k3d(Y,X,Z) ~= K_fluid && k3d(Y,X+1,Z) ~= K_fluid && k3d(Y,X-1,Z) == K_fluid && k3d(Y+1,X,Z) == K_fluid && k3d(Y-1,X,Z) ~= K_fluid && k3d(Y,X,Z+1) ~= K_fluid && k3d(Y,X,Z-1) == K_fluid)
                    T1(Y,X,Z) = (1/((1/(dX^2)+1/(dY^2)+1/(dZ^2))*k3d(Y,X,Z)+hc*(1/dX+1/dY+1/dZ)))*(hc*(T3d(Y,X-1,Z)/dX+T3d(Y+1,X,Z)/dY+T3d(Y,X,Z-1)/dZ)+(T3d(Y,X+1,Z)/(dX^2)+T3d(Y-1,X,Z)/(dY^2)+T3d(Y,X,Z+1)/(dZ^2))*k3d(Y,X,Z));
                    Transfer(Y,X,Z) = 'X <=, Y v, Z -';
                    Q1(Y,X,Z) = -hc*dY*dZ*(T3d(Y,X-1,Z) - T3d(Y,X,Z)) -hc*dX*dZ*(T3d(Y+1,X,Z) - T3d(Y,X,Z)) -hc*dX*dY*(T3d(Y,X,Z-1) - T3d(Y,X,Z));

                end                         
            end
        end
    end
    T3d = T1;
    end
    for M = (Xn_ins+2):(Xn_ins+Xn_max+1)
        for N = 3:Zn_tot-2 
        Q_bottom(1,M,N) = (dY*dZ*(T3d(2,M+1,N) - T3d(2,M,N)))/(dX/(k3d(2,M,N))) + (dY*dZ*(T3d(2,M-1,N) - T3d(2,M,N)))/(dX/(k3d(2,M,N))) + (dX*dZ*(T_boundaryZ(2,M,N) - T3d(2,M,N)))/(0.5*dY/(k3d(2,M,N))) + (dX*dY*(T3d(2,M,N+1) - T3d(2,M,N)))/(dZ/(k3d(2,M,N))) + (dX*dY*(T3d(2,M,N-1) - T3d(2,M,N)))/(dZ/(k3d(2,M,N)));
        end
    end

    Tc_1 = T_PV - 273;
    V = Voc +N_cells*kb*Tcel_ref*(log(Irr/I_ref));
    I = Isc*(Irr/I_ref);
    eta_E1 = ((I*V*FF)/(Irr*L*W))*(1+k_nu*(Tc_1-Tcel_ref));
    Quele1 = eta_E1*(Irr*(1-eta_glassc)*L*W); 

    Q_loss_bottom = abs(0.2*sum(Q_bottom(:)));         % sum(Q_air(:));
    Quth1 = 0.8*sum(Q_bottom(:));
    eta_tha1 = Quth1/(Irr*(1-eta_glassc)*L*W);

    Q_loss_top = Loss_convection + Loss_radiation;
    Q_in = Irr*(1-eta_glassc)*L*W;
    Q_out = Quth1 + Q_loss_top + Quele1 + Q_loss_bottom;
    N_correction = Q_in/Q_out; 

    Two_1 = Quth1/(C_w*dm_w)+Twin-273;
    if N_correction < 1.01 && N_correction > 0.99
        break
    end
    if N_correction < 0.5
        T_PV = T_PV - 1;
    elseif N_correction < 0.8 && N_correction > 0.5 
        T_PV = T_PV - 0.3;
    elseif N_correction < 0.9 && N_correction > 0.8 
        T_PV = T_PV - 0.1;
    elseif N_correction < 1 && N_correction > 0.9 
        T_PV = T_PV - 0.05;
    elseif N_correction > 1.5
        T_PV = T_PV + 1;
    elseif N_correction > 1.2 && N_correction < 1.5
        T_PV = T_PV + 0.3;    
    elseif N_correction > 1.1 && N_correction < 1.2
        T_PV = T_PV + 0.1;
    elseif N_correction > 1 && N_correction < 1.1
        T_PV = T_PV + 0.05;
    end
end

    % T2 = T3d - 273;
    % t3middle = T2(:,:,2);
    % Qmiddle = Transfer(:,:,2);
    % Q1middle = Q1(:,:,2);
    % outputs
    Two(1,hour) = Two_1; 
    Tc(1,hour) = Tc_1; 
    eta_E(1,hour) = eta_E1;
    eta_tha(hour) = eta_tha1;
    eta_tot(hour) = eta_E1 + eta_tha1;
    Quele(hour) = Quele1;
    Quth(hour) = Quth1;
    Power_loss_day(hour) = Q_loss_top;
    end
    Start_temp = Tc_1;
end

%% figures
% Coordinates
% Centermatrix_X = repmat((0:dX:Xn_tot*dX-dX),Yn_tot,1,Zn_tot); %X coordinate of the numeric centerpoints
% Centermatrix_Y = repmat((0:dY:Yn_tot*dY-dY)',1,Xn_tot,Zn_tot); %Y coordinate of the numeric centerpoints
% Centermatrix_Z = repmat(permute((0:dZ:Zn_tot*dZ-dZ),[1 3 2]),Yn_tot,Xn_tot,1); %Y coordinate of the numeric centerpoints
% 
% scrsz = get(0,'ScreenSize');
% figure('Position',[(scrsz(3)-figurew)/2,(scrsz(4)-figureh)/2,figurew, figureh]);
% xpl = Centermatrix_X(:);
% ypl = Centermatrix_Y(:);
% zpl = Centermatrix_Z(:);
% 
% minTemperature = min([T_ambient - 273,T_inflow - 273,T_PV - 273]);
% maxTemperature = max([T_ambient - 273,T_inflow - 273,T_PV - 273]);
% temperature = T2(:);
% numPoints = length(xpl);
% cmap = turbo(numPoints);
% percentage = (temperature-minTemperature) / (maxTemperature - minTemperature);
% indexes = round(percentage * (numPoints - 1)  + 1);
% tempColors = cmap(indexes, :);
% scatter3(xpl, ypl, zpl, 80, tempColors, 'filled');
% colormap(turbo)
% colorbar;
% caxis([minTemperature maxTemperature]);
% daspect([1 1 1]);
% 
% %middle slab
% temperature2 = reshape(T2(:,:,(Zn_tot)/2),[(Xn_tot)*(Yn_tot),1]);
% xpl2 = reshape(Centermatrix_X(:,:,(Zn_tot)/2),[(Xn_tot)*(Yn_tot),1]); 
% ypl2 = reshape(Centermatrix_Y(:,:,(Zn_tot)/2),[(Xn_tot)*(Yn_tot),1]); 
% numPoints2 = length(temperature2);
% cmap2 = turbo(numPoints2);
% percentage2 = (temperature2-minTemperature) / (maxTemperature - minTemperature);
% indexes2 = round(percentage2 * (numPoints2 - 1)  + 1);
% tempColors2 = cmap2(indexes2, :);
% 
% figure('Position',[(scrsz(3)-figurew)/2,(scrsz(4)-figureh)/2,figurew, figureh]); 
% scatter(xpl2,ypl2, 80, tempColors2, 'filled')
% colormap(turbo)
% colorbar;
% caxis([minTemperature maxTemperature]);
% daspect([4 1 1]);

%%
eta_E(Ia<=1) = 0;
eta_tha(eta_tha1<0) = 0;
eta_tha(Ia<=1) = 0;

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

pvt_collector_output.eta_tha=eta_tha1;
pvt_collector_output.eta_E=eta_E1;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.total_thermal_output=total_thermal_output;
pvt_collector_output.total_electrical_output=total_electrical_output;

%% Plot
period = WEATHER_output.Period;
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
ylabel('Temperature [^oC]');
yyaxis right
plot(xx,Ia,'DisplayName','I_{sun}')
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

f = figure(101);
plot(xx,eta_tha*100,'DisplayName','\eta_{th}')
hold on
plot(xx,eta_E*100,'DisplayName','\eta_{el}')
ylabel('Efficiency [%]');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Efficiency Plot');
    xlim([1 time])
    legend()
    SetFigureDefaults(f)

f = figure(102);
plot(xx,Quele,'DisplayName','P_{el}')
hold on
plot(xx,Quth,'DisplayName','P_{th}')
ylabel('Power [W]');
if time/24 >=31
        xticks(ticks); xticklabels(month_labelss)
    else
        xlabel('Time [hours]')
end
    title('Power Plot');
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
