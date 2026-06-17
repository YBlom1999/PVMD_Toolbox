function [V_module,I,parameters_1, parameters_2,parameters_3,indDiode] = TripleTandemModule(TOOLBOX_input,numCells, Acell,T_,J_,Irr,SubMod_ind)
% TripleTandemModule calculates the electrical performance of a triple tandem module
% and can simulate different module interconnections
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% numCells : double
%   The number of cells
% Acell: double
%   The area of the cell
% T_: double
%   The temperature of all cells
% J_: double
%   The photoreceived of all cells
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
%   The parameters of the equivalent circuit of the middle cell
% parameters_3 : double
%   The parameters of the equivalent circuit of the middle cell
% indDiode : double
%   The index of which bypass diodes are active
%
% Written by Y. Blom

%% Determine the necessary current range for the IV curves

N_submodules = length(Acell);
Imax = 0;
for Submod_i = 1:N_submodules
    Imax_new = max(J_{Submod_i},[],'all')*Acell(Submod_i)*1.05;
    Imax = max(Imax,Imax_new);
end
%Resolution of the IV curve in A
I=-2:5e-3:Imax;




%% Calculate cell IVs
Mod_ind = SubMod_ind(1); %Index of submodule
Cell_ind = sum(SubMod_ind(1:1)==Mod_ind); %Index of cell in submodule
Shading = TOOLBOX_input.electric.shading(Mod_ind);
IVtype = TOOLBOX_input.electric.IVtypeTop;
[V_1,day1,night1,parameters_1] = Simulated_IV(TOOLBOX_input,Acell(Mod_ind),numCells(Mod_ind),Shading,IVtype,J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind},Irr{Mod_ind},...
    I, 1,TOOLBOX_input.electric.degtype(1));

Mod_ind = SubMod_ind(2); %Index of submodule
Cell_ind = sum(SubMod_ind(1:2)==Mod_ind); %Index of cell in submodule
Shading = TOOLBOX_input.electric.shading(Mod_ind);
IVtype = TOOLBOX_input.electric.IVtypeMiddle;
[V_2,day2,night2,parameters_2] = Simulated_IV(TOOLBOX_input,Acell(Mod_ind),numCells(Mod_ind),Shading,IVtype,J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind},Irr{Mod_ind},...
    I, 2,TOOLBOX_input.electric.degtype(2));

Mod_ind = SubMod_ind(3); %Index of submodule
Cell_ind = sum(SubMod_ind(1:3)==Mod_ind); %Index of cell in submodule
Shading = TOOLBOX_input.electric.shading(Mod_ind);
IVtype = TOOLBOX_input.electric.IVtypeBot;
[V_3,day3,night3,parameters_3] = Simulated_IV(TOOLBOX_input,Acell(Mod_ind),numCells(Mod_ind),Shading,IVtype,J_{Mod_ind}(:,:,Cell_ind),T_{Mod_ind},Irr{Mod_ind},...
    I, 3,TOOLBOX_input.electric.degtype(3));

%% Compare day and night for both cells and give warning if necessary
if ~isequal(day1,day2,day3) || ~isequal(night1,night2,night3)
    warning('The calculated days and nights are not the same') 
end
day = day1;
night = night1;

%% Calculate modules IV

tripleType = TOOLBOX_input.electric.tripleType;
R_int = TOOLBOX_input.electric.resistance;
effLC = TOOLBOX_input.electric.LC_eff;
Nby = TOOLBOX_input.electric.numBypassDiodes;

if strcmp(tripleType,'ABC') %All series connected modules
    V_2 = includeLumCouplingTandem(V_1,nan,V_2,I,effLC(1));
    V_3 = includeLumCouplingTandem(V_1,V_2,V_3,I,effLC(2:3));

    V_tan = V_1+V_2+V_3 - R_int*I;
    if strcmp(TOOLBOX_input.electric.ModuleType,'Series')
        [V_module,I,indDiode] = connectCells2ModuleSeries(V_tan,Nby,numCells,I,day,night);
    elseif strcmp(TOOLBOX_input.electric.ModuleType,'ButterFly')
        [V_module,I,indDiode] = connectCells2ModuleButterfly(V_tan,Nby,numCells,I,day,night);
    end
elseif strcmp(tripleType,'A-BC') %Top cell seperate, middle and bottom cell together
    V_1 = V_1 - R_int(1)*I;
    if strcmp(TOOLBOX_input.electric.ModuleType{1},'Series')
        [V_module1,I,indDiode1] = connectCells2ModuleSeries(V_1,Nby(1),numCells(1),I,day1,night1);
    elseif strcmp(TOOLBOX_input.electric.ModuleType{1},'ButterFly')
        [V_module1,I,indDiode1] = connectCells2ModuleButterfly(V_1,Nby(1),numCells(1),I,day1,night1);
    end
    
    V_3 = includeLumCouplingTandem(V_2,nan,V_3,I,effLC(3));
    V_tan = V_2+V_3 - R_int(2)*I;
    if strcmp(TOOLBOX_input.electric.ModuleType{2},'Series')
        [V_module2,I,indDiode2] = connectCells2ModuleSeries(V_tan,Nby(2),numCells(2),I,day2,night2);
    elseif strcmp(TOOLBOX_input.electric.ModuleType{2},'ButterFly')
        [V_module2,I,indDiode2] = connectCells2ModuleButterfly(V_tan,Nby(2),numCells(2),I,day2,night2);
    end

    V_module=[V_module1;V_module2];
    indDiode=[indDiode1;indDiode2];
elseif strcmp(tripleType,'AB-C') %Top and middle cell together, bottom cell seperate
    V_2 = includeLumCouplingTandem(V_1,nan,V_2,I,effLC(1));
    V_tan = V_1+V_2 - R_int(1)*I;
    if strcmp(TOOLBOX_input.electric.ModuleType{1},'Series')
        [V_module1,I,indDiode1] = connectCells2ModuleSeries(V_tan,Nby(1),numCells(1),I,day1,night1);
    elseif strcmp(TOOLBOX_input.electric.ModuleType{2},'ButterFly')
        [V_module1,I,indDiode1] = connectCells2ModuleButterfly(V_tan,Nby(1),numCells(1),I,day1,night1);
    end

    V_3 = V_3 - R_int(2)*I;
    if strcmp(TOOLBOX_input.electric.ModuleType{2},'Series')
        [V_module2,I,indDiode2] = connectCells2ModuleSeries(V_3,Nby(2),numCells(2),I,day3,night3);
    elseif strcmp(TOOLBOX_input.electric.ModuleType{2},'ButterFly')
        [V_module2,I,indDiode2] = connectCells2ModuleButterfly(V_3,Nby(2),numCells(2),I,day3,night3);
    end

    V_module=[V_module1;V_module2];
    indDiode=[indDiode1;indDiode2];
end
end
