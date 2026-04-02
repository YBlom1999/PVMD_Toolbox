function [P_series,P_shunt,P_NRRI,P_NRRV,PowerRatio] = ElectricLossesTandem_avg(Elec_Losses_input,TOOLBOX_input,MODULE_output,CONSTANTS)
%ElectricLossesTandem_avg Calculates the electrical losses of the system
% with tandem modules.
%
% This function calculates the electrical losses for the PV module (tandem junction) on cell level. The
% losses are calculated as if the cell is operating on MPP
%
% Parameters
% ----------
% Elec_Losses_input : struc
%   The needed inputs for the electrical losses calculation
% TOOLBOX_input : struct
%   Simulation parameters
% MODULE_output : struct
%   Simulation results of the MODULE module
% CONSTANTS : struct
%   Physical constants
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
% PowerRatio : double
%   The power ratio of the top and bottom cell: Pmmp1/(Pmpp1+Pmpp2)
%
% Developed by Y. Blom


I_abs1 = Elec_Losses_input.I_abs1;
I_abs2 = Elec_Losses_input.I_abs2;
V_opt1 = Elec_Losses_input.V_opt1;
V_opt2 = Elec_Losses_input.V_opt2;
Parameters1 = Elec_Losses_input.Parameters1;
Parameters2 = Elec_Losses_input.Parameters2;
T_cell = Elec_Losses_input.T_cell;

% constants
q = CONSTANTS.q;
k = CONSTANTS.k;

N_cells = MODULE_output.N;
correction1 = I_abs1'-mean(Parameters1(:,:,1));
correction2 = I_abs2'-mean(Parameters2(:,:,1));

P_series1 = zeros(length(I_abs1),N_cells);
P_series2 = zeros(length(I_abs2),N_cells);
P_shunt1 = zeros(length(I_abs1),N_cells);
P_shunt2 = zeros(length(I_abs2),N_cells);
P_NRRV1 = zeros(length(I_abs1),N_cells);
P_NRRV2 = zeros(length(I_abs2),N_cells);
P_NRRI1 = zeros(length(I_abs1),N_cells);
P_NRRI2 = zeros(length(I_abs2),N_cells);

for i = 1:N_cells

    %Electrical properties
    Iph1 = Parameters1(i,:,1)'+correction1';
    Rs1 = Parameters1(i,:,2)';
    Rsh1 = Parameters1(i,:,3)';
    n1 = Parameters1(i,:,4)';
    I01 = Parameters1(i,:,5)';

    Iph2 = Parameters2(i,:,1)'+correction2';
    Rs2 = Parameters2(i,:,2)';
    Rsh2 = Parameters2(i,:,3)';
    n2 = Parameters2(i,:,4)';
    I02 = Parameters2(i,:,5)';
    [I_mpp1, V_mpp1] =  MPPfinder(Iph1,Rs1,Rsh1,n1,I01,T_cell);
    [I_mpp2, V_mpp2] =  MPPfinder(Iph2,Rs2,Rsh2,n2,I02,T_cell);
    V_mpp1(isnan(I_mpp1)) = 0;
    I_mpp1(isnan(I_mpp1)) = 0;
    V_mpp2(isnan(I_mpp2)) = 0;
    I_mpp2(isnan(I_mpp2)) = 0;

    PowerRatio = (I_mpp1.*V_mpp1)'./(I_mpp1.*V_mpp1+I_mpp2.*V_mpp2)';
    PowerRatio(isnan(PowerRatio)) = 0;

    V_series1 = I_mpp1.*Rs1';
    V_series2 = I_mpp1.*Rs2';

    V_NRRV1 = V_opt1-V_mpp1-V_series1;
    V_NRRV2 = V_opt2-V_mpp2-V_series2;

    I_shunt1 = (V_mpp1+I_mpp1*Rs1)./Rsh1';
    I_shunt2 = (V_mpp2+I_mpp1*Rs2)./Rsh2';
    I_shunt1(isnan(I_shunt1)) = 0;
    I_shunt2(isnan(I_shunt2)) = 0;
    I_shunt1(isinf(I_shunt1)) = 0;
    I_shunt2(isinf(I_shunt2)) = 0;

    I_NRRI1 = (Iph1'-I_mpp1-I_shunt1);
    I_NRRI2 = (Iph2'-I_mpp2-I_shunt2);

    P_series1(:,i) = I_mpp1'.*V_series1';
    P_series2(:,i) = I_mpp2'.*V_series2';

    P_shunt1(:,i) = I_shunt1'.*(V_mpp1+V_series1)';
    P_shunt2(:,i) = I_shunt2'.*(V_mpp2+V_series2)';

    P_NRRV1(:,i) = V_NRRV1'.*Iph1;
    P_NRRV2(:,i) = V_NRRV2'.*Iph2;

    P_NRRI1(:,i) = I_NRRI1'.*(V_mpp1+V_series1)';
    P_NRRI2(:,i) = I_NRRI2'.*(V_mpp2+V_series2)';


end
P_series = sum(P_series1+P_seri