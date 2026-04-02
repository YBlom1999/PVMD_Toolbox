function [P_con,P_mismatch,P_cable,P_inv] = SystemLossesTandem2T(Sys_Losses_input,TOOLBOX_input,ELECTRIC_output,MODULE_output)
%SystemLossesTandem Calculates the system losses of the PV system
% with tandem modules.
%
% This function calculates the system losses of the PV module for tandem
% cells
%
% Parameters
% ----------
% Sys_Losses_input : struc
%   The needed inputs for the system losses calculation
% TOOLBOX_input : struct
%   Simulation parameters
% ELECTRIC_output : struct
%   Simulation results of the ELECTRIC module
% MODULE_output : struct
%   Simulation results of the MODULE module
%
% Returns
% -------
% P_con : double
%   The losses due to the connection resistance
% P_mismatch : double
%   The losses due to current matching
% P_cable : double
%   The current loss due to cable losses
% P_inv : double
%   The voltage loss due to inverter losses
%
% Developed by Y. Blom

I_mpp = Sys_Losses_input.I_mpp;
V_mpp = Sys_Losses_input.V_mpp;
I_abs1 = Sys_Losses_input.I_abs1;
I_abs2 = Sys_Losses_input.I_abs2;
Parameters1 = Sys_Losses_input.Parameters1;
Parameters2 = Sys_Losses_input.Parameters2;
T_cell = Sys_Losses_input.T_cell;

N_cells = MODULE_output.N;
correction1 = I_abs1'-mean(Parameters1(:,:,1));
correction2 = I_abs2'-mean(Parameters2(:,:,1));


%Electrical properties
R_con = TOOLBOX_input.electric.resistance;
V_mod = V_mpp+I_mpp*R_con*N_cells;
P_con = N_cells*R_con*I_mpp.^2;

P_mpp_cell = zeros(length(V_mpp),N_cells);
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

    P_mpp_cell(:,i) = I_mpp1'.*V_mpp1'+I_mpp2'.*V_mpp2';

end
P_mismatch = sum(P_mpp_cell,2)-I_mpp.*V_mod;
if (isfield(TOOLBOX_input, 'runACConversionPart')==1 && TOOLBOX_input.runACConversionPart == 1)
    Pdc = Sys_Losses_input.Pdc;
    Pac = Sys_Losses_input.Pac;
    panels = TOOLBOX_input.Conversion.Parallel_Modules*TOOLBOX_input.Conversion.Series_Modules;
    if exist('TOOLBOX_input.Conversion.CableLossesDetailedCalculation','var') && (isnumeric(TOOLBOX_input.Conversion.CableLossesDetailedCalculation) == 1)
        P_cable = TOOLBOX_input.Conversion.CableLoss.Cable_loss_total;
    else
        P_cable = panels*Pdc*TOOLBOX_input.Conversion.CableLoss/100;
    end
    P_