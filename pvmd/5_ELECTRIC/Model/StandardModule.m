function [V_module,I,param,reconfig_setting,indDiode] = StandardModule(TOOLBOX_input,numCells, Acell,Temp,Jabs,Jabs_STC,Irr)
% StandardModule calculates the electrical performance of a standard module
% and can simulate series-connected, butterfly, and reconfigurable modules
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
% Jabs_STC : double
%   The absorbed current at STC
% Irr : double
%   The received irradiance of all cells
%
% Returns
% -------
% V_module : double
%   The voltage of the module IV curve
% I : double
%   The current of the module IV curve
% param : double
%   The parameters of the equivalent circuit
% reconfig_setting : struct
%   The settings of the reconfigurable modules
% indDiode : double
%   The index of which bypass diodes are active
%
% Written by M. R. Vogt
% Adjusted by Y. Blom


%maximum simulation current
Imax=max(Jabs,[],'all')*Acell*1.05;
%Resolution of the IV curve in A
I=0:(Imax/500):Imax;

%% Calculate cell IVs
if TOOLBOX_input.electric.InterpolateIV
    [Vcell,day,night] = InterpolateCell(Acell,Irr,Temp,TOOLBOX_input,numCells,I,Jabs_STC);
    param = nan;
else
    Shading = TOOLBOX_input.electric.shading;
    IVtype = TOOLBOX_input.electric.IVtype;
    [Vcell,day,night,param] = Simulated_IV(TOOLBOX_input, Acell, numCells, Shading, IVtype, Jabs, Temp, Irr, I,1,1);
end


%% Account for series resistance
R = TOOLBOX_input.electric.resistance;
Vcell=Vcell - I*R;

%% Calculate modules IV
if strcmp(TOOLBOX_input.electric.ModuleType,'Series')
    Nby = TOOLBOX_input.electric.numBypassDiodes;
    [V_module,I,indDiode] = connectCells2ModuleSeries(Vcell,Nby,numCells,I,day,night);
    reconfig_setting = nan;
elseif strcmp(TOOLBOX_input.electric.ModuleType,'Butterfly')
    Nby = TOOLBOX_input.electric.numBypassDiodes;
    [V_module,I,indDiode] = connectCells2ModuleButterfly(TOOLBOX_input,Vcell,Nby,numCells,I,day,night);
    reconfig_setting = nan;
elseif strcmp(TOOLBOX_input.electric.ModuleType,'REC')
    Rsw = TOOLBOX_input.electric.reconfig.Rsw;
    algo = TOOLBOX_input.electric.reconfig.algo;
    nc_r = TOOLBOX_input.electric.reconfig.nc_r;
    [V_module,reconfig_setting] =combineToModuleIV_ReconfigMod(Vcell,numCells,I,day,night,Rsw,algo,nc_r);
    indDiode = nan;
end
end