function [k_pred,factor] = DegMoistureIngress(WEATHER_output,THERMAL_output,ELECTRIC_output,TOOLBOX_input,T,Ea,A,C,Time,CONSTANTS)
%DegMoistureIngress Calculates the degradation due to moisture ingress
%
% This function calculates the degradation rate caused by the moisture
% ingress into the PV module
%
% Parameters
% ----------
% WEATHER_output: struct
%   Simulation result of the electrical part
% THERMAL_output: struct
%   Simulation result of the electrical part
% ELECTRIC_output: struct
%   Simulation result of the electrical part
% TOOLBOX_input : struct
%   Simulation parameters
% T_repeat : double
%   Temperature of the module
% Ea : double
%   The activation energy of discoloration
% A : double
%   The pre-exponential constant of discoloration
% C : double
%   The constant indicating how fast the module degradates.
% Time: double
%   The time for which the degradation needs to be simulated
% CONSTANTS : struct
%   Structure of physical constants
%
% Returns
% -------
% k_pred: double
%   The predicted degradation rate due to discoloration
%
% Developed by by Youri Blom
k_b = CONSTANTS.k_b;
q = CONSTANTS.q;
T_STC = CONSTANTS.T_STC;

%Determine RMC
if TOOLBOX_input.Degradation.loadCOMSOL
    [~,~,data_folder] = get_folder_structure;
    locations_dir = fullfile(data_folder,'Degradation','RMC');
    RMC_file = fullfile(locations_dir,TOOLBOX_input.Degradation.COMSOL_file);
    load(RMC_file,'RMC_data');
    Time_sim = RMC_data(:,1);
    RMC_sim = RMC_data(:,2);

    

    %The RMC needs to be repeated for the operating after 20 years
    if Time(end) > Time_sim(end)
        N_repeat = (Time(end)-Time_sim(end))/8760;
        t_1 = find(Time_sim==Time_sim(end)-8760)+1;
        t_2 = find(Time_sim==Time_sim(end));
        delta_t_sim = Time_sim(2)-Time_sim(1);
        RMC_repeat = RMC_sim(t_1:t_2);

        Time_sim = [Time_sim; (Time_sim(end)+delta_t_sim:delta_t_sim:Time(end))'];
        RMC_sim = [RMC_sim;repmat(RMC_repeat,N_repeat,1)];

    end

    RMC = max(interp1(Time_sim,RMC_sim,Time),0);
else


end

%Different equations can be used:
%Peck model
k_pred = A*exp(-q*Ea/k_b./T(Time)).*RMC.^C;
%Eyring model
% k_pred = A*exp(-q*Ea/k_b/T(Time)-C./RMC);
%Exponential model
% k_pred = A*exp(-q*Ea/k_b./T(Time)).*exp(C*RMC);
P_loss = sum(k_pred);

%Load Data file
filename = TOOLBOX_input.Degradation.MoistureDeg_file;
[~,~,data_folder] = get_folder_structure;
filename_path =fullfile(data_folder,'Degradation','MoistureIngress',filename);
load(filename_path,'Dosage','trend_I0','trend_Iph','trend_Rs','trend_Rsh','trend_n');

Iph_orig = ELECTRIC_output.Parameters_STC_1(1,1);
Rs_orig = ELECTRIC_output.Parameters_STC_1(1,2);
Rsh_orig = ELECTRIC_output.Parameters_STC_1(1,3);
n_orig = ELECTRIC_output.Parameters_STC_1(1,4);
I0_orig = ELECTRIC_output.Parameters_STC_1(1,5);

%Make array of power loss for different dosage to find which dossage is
%needed
P_range = zeros(1,length(Dosage));
Voltage = 0:0.01:1;
Vth=k_b*T_STC/q;
for i = 1:length(Dosage)
    Iph_new = Iph_orig*trend_Iph(i)/trend_Iph(1);
    I0_new = 10^(log10(I0_orig)*trend_I0(i)/trend_I0(1));
    Rs_new = Rs_orig*trend_Rs(i)/trend_Rs(1);
    Rsh_new = Rsh_orig*trend_Rsh(i)/trend_Rsh(1);
    n_new = n_orig*trend_n(i)/trend_n(1);
    
    z=(Rs_new*I0_new/(n_new*Vth*(1+Rs_new/Rsh_new)))*exp((Rs_new*(Iph_new+I0_new)+Voltage)./(n_new*Vth*(1+Rs_new/Rsh_new)));
    I_new=(Iph_new+I0_new-Voltage/(Rsh_new))/(1+Rs_new/Rsh_new)-lambertw(z).*(n_new*Vth)/Rs_new;

    P_range(i) = max(I_new.*Voltage);
    if i>1 && (P_range(i) > P_range(i-1))
        P_range(i) = P_range(i-1)-0.0001;
    end
end
P_loss_range = 1-P_range/P_range(1);
Dosage_needed = interp1(P_loss_range,Dosage,P_loss);

%Update the parameters
factor_Iph = interp1(Dosage,trend_Iph,Dosage_needed)/interp1(Dosage,trend_Iph,0);
factor_I0 = interp1(Dosage,trend_I0,Dosage_needed)/interp1(Dosage,trend_I0,0);
factor_Rs = interp1(Dosage,trend_Rs,Dosage_needed)/interp1(Dosage,trend_Rs,0);
factor_Rsh = interp1(Dosage,trend_Rsh,Dosage_needed)/interp1(Dosage,trend_Rsh,0);
factor_n = interp1(Dosage,trend_n,Dosage_needed)/interp1(Dosage,trend_n,0);

factor = [factor_Iph,factor_Rs,factor_Rsh,factor_n,factor_I0];
end


