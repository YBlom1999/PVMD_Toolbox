function [factors] = FactorsSilicon(Parameters,k_needed,T,filename)
%FactorsSilicon file for the degradation module of tandem modules
%
% This function calculates the degradation rate of a tandem module
% 
%
% Parameters
% ----------
% Parameters : double
%   Parameters of the equivalent circuit
% k_needed : double
%   The degradation rate that is needed.
% T : double
%   The temperature of the module.
% filename : double
%   The name of the moisture ingress degraation file
%
% Returns
% -------
% factors : double
%   The factors by which the parameters should be updated
%
% Developed by by Youri Blom

K=1.3806e-23;
q=1.6022e-19;
Vth=K*T/q;

[~,~,data_folder] = get_folder_structure;
filename_path =fullfile(data_folder,'Degradation','MoistureIngress',filename);
load(filename_path,'Dosage','trend_I0','trend_Iph','trend_Rs','trend_Rsh','trend_n');

Iph_orig = Parameters(1,1);
Rs_orig = Parameters(1,2);
Rsh_orig = Parameters(1,3);
n_orig = Parameters(1,4);
I0_orig = Parameters(1,5);

%Make array of power loss for different dosage to find which dossage is
%needed
P_range = zeros(1,length(Dosage));
Voltage = 0:0.01:1;
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
[~,un_ind] = unique(P_loss_range);
Dosage_needed = interp1(P_loss_range(un_ind),Dosage(un_ind),k_needed,'linear','extrap');

%Update the parameters
factor_Iph = interp1(Dosage,trend_Iph,Dosage_needed,'linear','extrap')/interp1(Dosage,trend_Iph,0);
factor_I0 = interp1(Dosage,trend_I0,Dosage_needed,'linear','extrap')/interp1(Dosage,trend_I0,0);
factor_Rs = interp1(Dosage,trend_Rs,Dosage_needed,'linear','extrap')/interp1(Dosage,trend_Rs,0);
factor_Rsh = interp1(Dosage,trend_Rsh,Dosage_needed,'linear','extrap')/interp1(Dosage,trend_Rsh,0);
factor_n = interp1(Dosage,trend_n,Dosage_needed,'linear','extrap')/interp1(Dosage,trend_n,0);

factors = [factor_Iph,factor_Rs,factor_Rsh,factor_n,factor_I0];

end