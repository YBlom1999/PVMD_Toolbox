function [CELL_FRONT,Lay,Int,absmat] = simulateAbsorptionFront(Lay,Int,absmat,SET,wav,absorptanceAll)
% simulateAbsorptionFront simulates the front side of the module
%
% This function calculates the absorption profile of the front side
%
% Parameters
% ----------
% Lay: struct
%   The layers of the simulated cell
% Int: struct
%   The interfaces between the layers
% Absmat: double
%   Index of which layers or interfaces are absorber materials
% rback: double
%   Reflectance of backsheet (nan for bifacial modules)
% Type: string
%   The technology type that is used for the module
% Submod_ind: double
%   Index of which absorber materials belong to which submodule
% SET: struct
%   An overview of the settings that are used
% wav: double
%   The wavelengths that are simulated
% absorptanceAll: boolean
%   Indicator of the absorptance of all layers should be used as output
%
% Returns
% -------
% CELL_FRONT : struct
%   The output of the front side of the cell
% Lay: struct
%   The layers of the simulated cell
% Int: struct
%   The interfaces between the layers
% Absmat: double
%   Index of which layers or interfaces are absorber materials

%
% Developed by Youri Blom

[Lay,Int,out] = GENPRO_core(Lay,Int,SET);    %run Genpro
nonTextLayers = find(~startsWith(out.Layer,'*'));
absmat = [1,nonTextLayers(absmat)'+1,length(out.Layer)];

%if the absorptance for all materials is needed, ix is changed
if absorptanceAll
    AbsorberMaterials = absmat;
    absmat = 1:size(out,1);
end
Jph_STC = zeros(SET.nr_ang_interval,length(absmat)-2);
Jph_STC(1,:) = table2array(out(absmat(2:end-1),4));      %[mA/cm2]

RAT_f = zeros(length(wav),SET.nr_ang_interval,length(absmat));
RAT_f(:,1,:) = table2array(out(absmat,5:end))'; %save output (first (R), absorber layers (A) and final (T))


%---INCIDENCE FROM THE FRONT---
for c = 2:SET.nr_ang_interval       %for every other angular interval on front side
    SET.incident_ang_interval = c;        %set incident angular interval
    [~,~,out] = GENPRO_core(Lay,Int,SET); %run
    RAT_f(:,c,:) = table2array(out(absmat,5:end))';           %put absorptance output in array
    Jph_STC(c,:) = table2array(out(absmat(2:end-1),4));       %[mA/cm2]
end
%---put output in struct and copy to workspace---
CELL_FRONT.RAT = RAT_f;
CELL_FRONT.wav = wav;
CELL_FRONT.aoi = ((1:SET.nr_ang_interval)'-0.5)*90/SET.nr_ang_interval;
CELL_FRONT.lay = out.Layer(absmat);
CELL_FRONT.Jph = Jph_STC;
%If all absorptance is needed, the index of the absorber materials is
%needed.
if absorptanceAll
    CELL_FRONT.Absmat_ind = AbsorberMaterials;
end

end