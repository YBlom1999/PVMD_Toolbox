function [CELL_output] = CELL_main(TOOLBOX_input)
% CELL_main Main file for the CELL module in the PVMD toolbox
%
% This function calculates the absorption profile of a PV module
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
%
% Returns
% -------
% CELL_output : struct
%   Output of the CELL block
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by unknown, editted by Youri Blom

if TOOLBOX_input.deviceOptic.runGenPro==false
    load(TOOLBOX_input.deviceOptic.GenProFile,'CELL_output');
    if isfield(CELL_output.CELL_FRONT,'Absmat_ind')
        TOOLBOX_input.deviceOptic.exportAbsorptanceAll = 1;
    else
        TOOLBOX_input.deviceOptic.exportAbsorptanceAll = 0;
    end
else
    [Lay, Int, absmat, rback, Type, Submod_ind, SET, wav] = prepareInputCell(TOOLBOX_input);
    
    disp('CELL calculation started. This can take a few minutes...')
    absorptanceAll = TOOLBOX_input.deviceOptic.exportAbsorptanceAll;
    [CELL_FRONT,Lay,Int,absmat] = simulateAbsorptionFront(Lay,Int,absmat,SET,wav,absorptanceAll);
    
    %---INCIDENCE FROM THE REAR---
    if isnan(rback)               %if rback is empty, bificial module is used
        CELL_REAR = simulateAbsorptionRear(Lay,Int,absmat,SET,wav,absorptanceAll);
    else                            %otherwise use fixed reflectance value rback
        CELL_REAR = [rback,0];
    end
    
    CELL_output.CELL_FRONT = CELL_FRONT;    %put output in single struct
    CELL_output.CELL_REAR = CELL_REAR;
    CELL_output.TYPE=Type;
    CELL_output.SUBMOD_IND = Submod_ind;
    %---
    disp('Wavelength and angle dependent absoptance calculated.')
end
end
