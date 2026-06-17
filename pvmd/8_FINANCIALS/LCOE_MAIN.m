function [LCOE_output] = LCOE_MAIN(TOOLBOX_input,ELECTRIC_output,DEGRADATION_output,CONVERSION_output)
% LCOE_MAIN calculates the costs of the PV system.
%
% This function calculates the LCOE or the NPC costs of the PV system. For
% residential PV systems it is possible to include heat pump and battery
% storage.
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% ELECTRIC_output : struct
%   Output of the ELECTRIC block
%
% Returns
% -------
% LCOE_output : struct
%   Output of this module block
%
% Developed by ?, editted by Youri Blom

if~TOOLBOX_input.script
    Str={'Select the situation'} ;
    S={'PV plant','Residential PV'};
    Situation=listdlg('PromptString',Str,'ListString',S,'SelectionMode','single','ListSize',[300 70]);
    TOOLBOX_input.FinancialPart.Situation = Situation;
end

if TOOLBOX_input.FinancialPart.Situation == 1
    [LCOE_output, TOOLBOX_input] = LCOE_PVplant(TOOLBOX_input,ELECTRIC_output,DEGRADATION_output);
elseif TOOLBOX_input.FinancialPart.Situation == 2
    [LCOE_output, TOOLBOX_input] = LCOE_residentialPV(TOOLBOX_input,ELECTRIC_output,DEGRADATION_output,CONVERSION_output);
end

disp('LCOE calculation finished.')
end