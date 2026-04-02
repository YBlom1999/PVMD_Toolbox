function [V_module,I,param,reconfig_setting,ind_Diode] = StandardModule(TOOLBOX_input,numCells, Acell,T_,J_,Jph_STC,Irr,Losses)
%%cordiantes the the tranistion from Temperatures, Generated current to
%over cell IV curves to module IV curves
%Author: Malte Vogt

%% Account for Isc given in data sheet model
if TOOLBOX_input.electric.runDatasheet
    JphGenPro=Jph_STC;
    if strcmp(TOOLBOX_input.electric.TYPE,'Butterfly')
        Conv_Datasheet_GenPro=TOOLBOX_input.electric.datasheetValues(2)./(JphGenPro(1)*Acell*2);
    else
        Conv_Datasheet_GenPro=TOOLBOX_input.electric.datasheetValues(2)./(JphGenPro(1)*Acell);
    end
else
    Conv_Datasheet_GenPro=1;
end


%maximum simulation current
Imax=max(J_,[],'all')*Acell*Conv_Datasheet_GenPro*1.05;
%Resolution of the IV curve in A
if isfield(TOOLBOX_input.electric,'TaylorParam')
%     I=0:(Imax/500):TOOLBOX_input.electric.TaylorParam.I0max/2;
    I = linspace(0,TOOLBOX_input.electric.TaylorParam.I0max/2,501);
else
    I=0:(Imax/500):Imax;
end

if strcmp(TOOLBOX_input.electric.TYPE,'Butterfly')
    ButterFly = 1;
else
    ButterFly = 0;
end
%% Calculate cell IVs
if TOOLBOX_input.electric.runDatasheet
    [V_,day,night,param] = Module_Datasheet(Acell,J_,T_,TOOLBOX_input.electric.shading,...
        numCells,I, TOOLBOX_input.electric.datasheetValues,TOOLBOX_input.electric,ButterFly,JphGenPro(1));
elseif TOOLBOX_input.electric.InterpolateIV
    [V_,day,night] = InterpolateCell(Acell,Irr,T_,TOOLBOX_input,numCells,I,Jph_STC);
    param = nan;
else
    [V_,day,night,param] = Simulated_IV(Acell,J_,T_, TOOLBOX_input.electric.shading,...
        numCells,I, TOOLBOX_input.electric.IVtype,TOOLBOX_input.electric,1,Losses); %Changed by youri
end


%% Calculate modules IV
R = TOOLBOX_input.electric.resistance;
if strcmp(TOOLBOX_input.electric.TYPE,'Non-REC') || strcmp(TOOLBOX_input.electric.TYPE,'Butterfly')
    Nby = TOOLBOX_input.electric.numBypassDiodes;
    [V_module,I,ind_Diode] = combineToModuleIV_StandardMod(V_,Nby,numCells,I,day,night,R,ButterFly,TOOLBOX_input);
    reconfig_setting = nan;
elseif strcmp(TOOLBOX_input.electric.TYPE,'REC')
    Rsw = TOOLBOX_input.electric.reconfig.Rsw;
    algo = TOOLBOX_input.electric.reconfig.algo;
    nc_r = TOOLBOX_input.electric.reconfig.nc_r;
    [V_module,reconfig_setting] =combineToModuleIV_ReconfigMod(V_,numCells,I,day,night,R,Rsw,algo,nc_r);
    ind_Diode = nan;
end



end
