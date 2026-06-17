function [V_module,I,parameters_1, parameters_2,indDiode] = TandemModule(TOOLBOX_input,numCells, Acell,Temp,Jabs,Irr,SubMod_ind)
% TandemModule calculates the electrical performance of a tandem module
% and can simulate 2T, 3T, and 4T modules
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% numCells : double
%   The number of cells
% Acell: double
%   The area of the cell
% Temp: double
%   The temperature of all cells
% Jabs: double
%   The absorbed current in all cells
% Irr : double
%   The received irradiance of all cells
% Submod_ind : double
%   The index of in which submodules, the cells are located
%
% Returns
% -------
% V_module : double
%   The voltage of the module IV curve
% I : double
%   The current of the module IV curve
% parameters_1 : double
%   The parameters of the equivalent circuit of the top cell
% parameters_2 : double
%   The parameters of the equivalent circuit of the bottom cell
% indDiode : double
%   The index of which bypass diodes are active
%
% Written by M. R. Vogt
% Adjusted by Y. Blom

%% Determine the necessary current range for the IV curves

%maximum simulation current
N_submodules = length(Acell);
Imax = 0;
for Submod_i = 1:N_submodules
    Imax_new = max(Jabs{Submod_i},[],'all')*Acell(Submod_i)*1.05;
    Imax = max(Imax,Imax_new);
end
%Resolution of the IV curve in A
I=0:5e-3:Imax;



%% Calculate cell IVs
Mod_ind = SubMod_ind(1); %Index of submodule
Cell_ind = sum(SubMod_ind(1:1)==Mod_ind); %Index of cell in submodule
Shading = TOOLBOX_input.electric.shading(Mod_ind);
IVtype = TOOLBOX_input.electric.IVtypeTop;
[V_1,day1,night1,parameters_1] = Simulated_IV(TOOLBOX_input,Acell(Mod_ind),numCells(Mod_ind),Shading,IVtype,Jabs{Mod_ind}(:,:,Cell_ind),Temp{Mod_ind},Irr,...
    I, 1,TOOLBOX_input.electric.degtype(1));


Mod_ind = SubMod_ind(2); %Index of submodule
Cell_ind = sum(SubMod_ind(1:2)==Mod_ind); %Index of cell in submodule
Shading = TOOLBOX_input.electric.shading(Mod_ind);
IVtype = TOOLBOX_input.electric.IVtypeBot;
[V_2,day2,night2,parameters_2] = Simulated_IV(TOOLBOX_input,Acell(Mod_ind),numCells(Mod_ind),Shading,IVtype,Jabs{Mod_ind}(:,:,Cell_ind),Temp{Mod_ind},Irr,...
    I, 2,TOOLBOX_input.electric.degtype(2));



%% Compare day and night for both cells and give warning if necessary
if ~isequal(day1,day2) || ~isequal(night1,night2)
    warning('The calculated days and nights are not the same') 
end
day = day1;
night = night1;

%% Calculate modules IV
R_int = TOOLBOX_input.electric.resistance;
effLC = TOOLBOX_input.electric.LC_eff;
Nby = TOOLBOX_input.electric.numBypassDiodes;


if TOOLBOX_input.electric.Terminals==2
    V_2 = includeLumCouplingTandem(V_1,nan,V_2,I,effLC);
    V_tan = V_1+V_2-R_int*I;
    
    if strcmp(TOOLBOX_input.electric.ModuleType,'Series')
        [V_module,I,indDiode] = connectCells2ModuleSeries(V_tan,Nby,numCells,I,day,night);
    elseif strcmp(TOOLBOX_input.electric.ModuleType,'Butterfly')
        [V_module,I,indDiode] = connectCells2ModuleButterfly(TOOLBOX_input,V_tan,Nby,numCells,I,day,night);
    end
elseif TOOLBOX_input.electric.Terminals==3
    VM_ratio_m = TOOLBOX_input.electric.VM_ratio_m;
    VM_ratio_n = TOOLBOX_input.electric.VM_ratio_n;
    configuration = TOOLBOX_input.electric.configuration;

    V_1 = V_1 - R_int(1)*I;
    V_2 = V_2 - R_int(2)*I;
    [V_module,I] = combineToModuleIV_3Terminal_literature(V_1,V_2,Nby,...
    numCells,I,day,night,configuration,VM_ratio_m,VM_ratio_n);
elseif TOOLBOX_input.electric.Terminals==4
    V_1 = V_1 - R_int(1)*I;
    if strcmp(TOOLBOX_input.electric.ModuleType{1},'Series')
        [V_module1,I,indDiode1] = connectCells2ModuleSeries(V_1,Nby(1),numCells(1),I,day1,night1);
    elseif strcmp(TOOLBOX_input.electric.ModuleType{1},'ButterFly')
        [V_module1,I,indDiode1] = connectCells2ModuleButterfly(TOOLBOX_input,V_1,Nby(1),numCells(1),I,day1,night1);
    end
    V_2 = V_2 - R_int(2)*I;
    if strcmp(TOOLBOX_input.electric.ModuleType{2},'Series')
        [V_module2,I,indDiode2] = connectCells2ModuleSeries(V_2,Nby(2),numCells(2),I,day2,night2);
    elseif strcmp(TOOLBOX_input.electric.ModuleType{2},'ButterFly')
        [V_module2,I,indDiode2] = connectCells2ModuleButterfly(TOOLBOX_input,V_2,Nby(2),numCells(2),I,day2,night2);
    end

    V_module=[V_module1;V_module2];
    indDiode=[indDiode1;indDiode2];
end
end
