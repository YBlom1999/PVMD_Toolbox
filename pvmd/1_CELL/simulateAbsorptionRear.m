function CELL_REAR = simulateAbsorptionRear(Lay,Int,absmat,SET,wav,absorptanceAll)
% simulateAbsorptionRear simulates the rear side of the module
%
% This function calculates the absorption profile of the rear side
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
% CELL_REAR : struct
%   The output of the front side of the cell
%
% Developed by Youri Blom

RAT_r = zeros(length(wav),SET.nr_ang_interval,length(absmat));
Jph_STC_r = zeros(SET.nr_ang_interval,length(absmat)-2);
for c = 1:SET.nr_ang_interval             %for every angular interval on rear side
    SET.incident_ang_interval =  SET.nr_ang_interval * (4 * length(Int) - 2) + c;      %set incident angular interval (counting back from final one)
    [~,~,out] = GENPRO_core(Lay,Int,SET); %run
    RAT_r(:,c,:) = table2array(out(absmat,5:end))';           %put absorptance output in array
    Jph_STC_r(c,:) = table2array(out(absmat(2:end-1),4));       %[mA/cm2]
    %CELL_REAR.Siai_test(c)=S.iai;
end
%---put output in struct and copy to workspace---
CELL_REAR.RAT = RAT_r;
CELL_REAR.wav = wav;
CELL_REAR.aoi = ((1:SET.nr_ang_interval)'-0.5)*90/SET.nr_ang_interval;
CELL_REAR.lay = out.Layer(absmat);
CELL_REAR.Jph = Jph_STC_r;
if absorptanceAll
    CELL_REAR.Absmat_ind = AbsorberMaterials;
end

end