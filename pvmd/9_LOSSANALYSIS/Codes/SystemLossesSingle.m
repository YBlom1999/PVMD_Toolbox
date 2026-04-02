function [P_con,P_mismatch,P_cable,P_inv] = SystemLossesSingle(Sys_Losses_input,TOOLBOX_input,ELECTRIC_output,MODULE_output)
%SystemLossesSingle Calculates the system losses of the system 
% with single junction modules.
%
% This function calculates the system losses of the PV module for single 
% junction cells
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
I_abs = Sys_Losses_input.I_abs;
Parameters = Sys_Losses_input.Parameters;
T_cell = Sys_Losses_input.T_cell;


N_cells = MODULE_output.N;
correction = I_abs'-mean(Parameters(:,:,1));
%Electrical properties
R_con = TOOLBOX_input.electric.resistance;
V_mod = V_mpp+I_mpp*R_con*N_cells;
P_con = N_cells*R_con*I_mpp.^2;

P_mpp_cell = zeros(length(V_mpp),N_cells);
for i = 1:N_cells

    %Electrical properties
    Iph = Parameters(i,:,1)'+correction';
    Rs = Parameters(i,:,2)';
    Rsh = Parameters(i,:,3)';
    n = Parameters(i,:,4)';
    I0 = Parameters(i,:,5)';

    [I_mpp_cell, V_mpp_cell] =  MPPfinder(Iph,Rs,Rsh,n,I0,T_cell);
    V_mpp_cell(isnan(I_mpp_cell)) = 0;
    I_mpp_cell(isnan(I_mpp_cell)) = 0;
    P_mpp_cell(:,i) = I_mpp_cell'.*V_mpp_cell';

end
P_mismatch = sum(P_mpp_cell,2)-I_mpp.*V_mod;

if (isfield(TOOLBOX_input, 'runACConversionPart')==1)
    if TOOLBOX_input.runACConversionPart == 1
        Pdc = Sys_Losses_input.Pdc;
        Pac = Sys_Losses_input.Pac;
        panels = TOOLBOX_input.Conversion.Parallel_Modules*TOOLBOX_input.Conversion.Series_Modules;
        if TOOLBOX_input.Conversion.CableLossesDetailedCalculation == 0
            P_cable = panels*Pdc*TOOLBOX_input.Conversion.CableLoss/100;
        else
            P_cable = TOOLBOX_input.Conversion.CableLoss.Cable_loss_total;
        end
        P_inv = Pdc*panels - P_cable-Pac;
    else
        P_cable = 