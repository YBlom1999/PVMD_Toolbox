function [P_series,P_shunt,P_NRRI,P_NRRV] = ElectricLossesSingle_avg(Elec_Losses_input,MODULE_output)
%ElectricLossesSingle_avg Calculates the electrical losses of the system 
% with single junction modules.
%
% This function calculates the electrical losses for the PV module (single junction) on cell level. The
% losses are calculated as if the cell is operating on MPP
%
% Parameters
% ----------
% Elec_Losses_input : struc
%   The needed inputs for the electrical losses calculation
% MODULE_output : struct
%   Simulation results of the MODULE module
%
% Returns
% -------
% P_series : double
%   The losses due to the series resistance
% P_shunt : double
%   The losses due to the shunt resistance
% P_NRRI : double
%   The current loss due to non radiative recombination
% P_NRRV : double
%   The voltage loss due to non radiative recombination
%
% Developed by Y. Blom

I_abs = Elec_Losses_input.I_abs;
V_opt1 = Elec_Losses_input.V_opt1;
Parameters = Elec_Losses_input.Parameters;
T_cell = Elec_Losses_input.T_cell;

N_cells = MODULE_output.N;
correction = I_abs'-mean(Parameters(:,:,1));

P_series = zeros(length(I_abs),N_cells);
P_shunt = zeros(length(I_abs),N_cells);
P_NRRV = zeros(length(I_abs),N_cells);
P_NRRI = zeros(length(I_abs),N_cells);

for i = 1:N_cells

    %Electrical properties
    Iph = Parameters(i,:,1)'+correction';
    Rs = Parameters(i,:,2)';
    Rsh = Parameters(i,:,3)';
    n = Parameters(i,:,4)';
    I0 = Parameters(i,:,5)';

    [I_mpp, V_mpp] =  MPPfinder(Iph,Rs,Rsh,n,I0,T_cell);
    V_mpp(isnan(I_mpp)) = 0;
    I_mpp(isnan(I_mpp)) = 0;

    V_series = I_mpp.*Rs';

    V_NRRV = V_opt1-V_mpp-V_series;

    I_shunt = (V_mpp+I_mpp*Rs)./Rsh';
    I_shunt(isnan(I_shunt)) = 0;
    I_NRRI = (Iph'-I_mpp-I_shunt);

    P_series(:,i) = I_mpp'.*V_series';
    P_shunt(:,i) = I_shunt'.*(V_mpp+V_series)';
    P_NRRV(:,i) = V_NRRV'.*Iph;
    P_NRRI(:,i) = I_NRRI'.*(V_mpp+V_series)';


end
P_series = sum(P_series,2);
P_shunt 