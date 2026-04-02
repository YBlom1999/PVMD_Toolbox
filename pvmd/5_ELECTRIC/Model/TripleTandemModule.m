function [V_module,I,parameters_1, parameters_2,parameters_3] = TripleTandemModule(TOOLBOX_input,numCells, Acell,T_,J_,Jph_STC,SubMod_ind) %Changed by youri
%%cordiantes the the tranistion from Temperatures, Generated current to
%over cell IV curves to module IV curves
%Author: Malte Vogt

%% Determine the necessary current range for the IV curves

%Account for Isc given in data sheet model
if TOOLBOX_input.electric.runDatasheet
    JscGenPro=Jph_STC;
    Conv_top=TOOLBOX_input.electric.datasheetValuesTop(2)/(JscGenPro(1)*Acell); 
    Conv_bot=TOOLBOX_input.electric.datasheetValuesBot(2)/(JscGenPro(2)*Acell);
    Conv_Datasheet_GenPro=max(Conv_top,Conv_bot);
else
    Conv_Datasheet_GenPro=1;
end    

N_submodules = length(Acell);
Imax = 0;
for Submod_i = 1:N_submodules
    Imax_new = max(J_{Submod_i},[],'all')*Acell(Submod_i)*Conv_Datasheet_GenPro*1.05;
    Imax = max(Imax,Imax_new);
end
%Resolution of the IV curve in A
I=-2:5e-3:Imax;




%% Calculate cell IVs
Mod_ind = SubMod_ind(1); %Index of submodule
Cell_ind = sum(SubMod_ind(1:1)==Mod_ind); %Index of cell in submodule
if TOOLBOX_input.electric.runDatasheet
    [V_1,day,night,parameters_1] = Module_Datasheet(Acell(Mod_ind),J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind},TOOLBOX_input.electric.shading(Mod_ind),...
        numCells(Mod_ind),I, TOOLBOX_input.electric.datasheetValuesTop); %Changed by Youri
else
    [V_1,day,night,parameters_1] = Simulated_IV(Acell(Mod_ind),J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind}, TOOLBOX_input.electric.shading(Mod_ind),...
        numCells(Mod_ind),I, TOOLBOX_input.electric.IVtypeTop,ones(1,5)); %Changed by youri
end

Mod_ind = SubMod_ind(2); %Index of submodule
Cell_ind = sum(SubMod_ind(1:2)==Mod_ind); %Index of cell in submodule
if TOOLBOX_input.electric.runDatasheet
    [V_2,day,night,parameters_2] = Module_Datasheet(Acell(Mod_ind),J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind},TOOLBOX_input.electric.shading(Mod_ind),...
        numCells(Mod_ind),I, TOOLBOX_input.electric.datasheetValuesMid); %Changed by Youri
else
    [V_2,day,night,parameters_2] = Simulated_IV(Acell(Mod_ind),J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind}, TOOLBOX_input.electric.shading(Mod_ind),...
        numCells(Mod_ind),I, TOOLBOX_input.electric.IVtypeMid,ones(1,5)); %Changed by youri
end

Mod_ind = SubMod_ind(3); %Index of submodule
Cell_ind = sum(SubMod_ind(1:3)==Mod_ind); %Index of cell in submodule
if TOOLBOX_input.electric.runDatasheet
    [V_3,day,night,parameters_3] = Module_Datasheet(Acell(Mod_ind),J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind},TOOLBOX_input.electric.shading(Mod_ind),...
        numCells(Mod_ind),I, TOOLBOX_input.electric.datasheetValuesBot); %Changed by Youri
else
    [V_3,day,night,parameters_3] = Simulated_IV(Acell(Mod_ind),J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind}, TOOLBOX_input.electric.shading(Mod_ind),...
        numCells(Mod_ind),I, TOOLBOX_input.electric.IVtypeBot,ones(1,5)); %Changed by youri
end


%% To Do: Compare day and night for both cells and give warning if necessary



%% Calculate modules IV

R_int = TOOLBOX_input.electric.resistance;
LC_eff = TOOLBOX_input.electric.LC_eff;
Type = TOOLBOX_input.electric.Type;
N_by = TOOLBOX_input.electric.numBypassDiodes;
[V_module] = combineToModuleIV_TripleTandem(V_1,V_2,V_3,N_by,numCells,I,day,night,R_int,LC_eff,Type,SubMod_ind);
  

end
