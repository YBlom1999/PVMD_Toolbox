function [ELECTRIC_output, TOOLBOX_input] = ELECTRIC_main(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output)
%ELECTRIC_MAIN Main file for the Electrical module in the PVMD toolbox
%
% This function calculates the module IV curve based on the generated
% charge carries and the temperature
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Simulation results of the CELL module
% MODULE_output : struct
%   Simulation results of the MODULE module
% WEATHER_output : struct
%   Simulation results of the WEATHER module
% THERMAL_output : struct
%   Simulation results of the THERMAL module
%
% Returns
% -------
% THERMAL_output : struct
%   Simulation results of the electric module
% TOOLBOX_input : struct
%   Update of TOOLBOX_input by adding simulation parameters of electric module
%
% Developed by by Abdallah Nour El Din and Malte Vogt

%%This function calculates the module IV curve based on the generated
%charge carries and the temperature
%The code was developed by Abdallah Nour El Din and Malte Vogt



%load input values
addpath(genpath("Model"));
TOOLBOX_input.runElectricPart = true;


if TOOLBOX_input.script
    Rsw = nan;
    algo = nan;
    nc_r = nan;
end
%including warnings for large data-sets
days_of_simulation = WEATHER_output.Period;

start_day = days_of_simulation(1);
start_month = days_of_simulation(2);
end_day = days_of_simulation(3);
end_month = days_of_simulation(4);

date1 = datetime(2005,start_month,start_day); %starting day of simulations
date2 = datetime(2005,end_month,end_day); %ending day of simulations
dur = hours(date2 - date1); %duration in hours

if TOOLBOX_input.electric.IncludeDissapatingHeat
    N_iter = TOOLBOX_input.electric.N_iter;
else
    N_iter = 1;
end

%keeping 31 as benchmark to describe small simulations
%Any simulation longer than 1 month/31 days is considered to be large
if ~TOOLBOX_input.runPeriodic
    if dur > 744 %if dur > number of hours in a month with 31 days
        warning('Simulation time will be too large for Datasheet Values!!!');
        disp('Recommended to proceed with Model Cell Values for DC Electrical Simulations');
    end
end

TYPE=CELL_output.TYPE;
SUBMOD_IND = CELL_output.SUBMOD_IND;
N_submod = max(SUBMOD_IND);


if TOOLBOX_input.runPeriodic
    Acell=MODULE_output.A;
    numCells=MODULE_output.N;
    Irr = WEATHER_output.Irr;
    T_ = cell(1,N_submod);
    J_ = cell(1,N_submod);
    for SubMod_i = 1:N_submod
        T_{SubMod_i}=THERMAL_output.T{SubMod_i} +273.15;
        if TOOLBOX_input.deviceOptic.exportAbsorptanceAll
            ind = CELL_output.CELL_FRONT.Absmat_ind(2:end-1)-2;
            J_{SubMod_i}=squeeze(WEATHER_output.J{SubMod_i}(:,:,ind))*1.6022e-19;
        else
            J_{SubMod_i}=WEATHER_output.J{SubMod_i}*1.6022e-19;
        end
    end
else %non-periodic simulation
    N_panels = length(MODULE_output.Panels);
    Acell=MODULE_output.Panels(1).Acell;
    numCells=MODULE_output.Panels(1).Ncells;
    Irr = cell(1,N_panels);
    T_ = cell(1,N_panels);
    J_ = cell(1,N_panels);
    for panel_i = 1:N_panels
        T_{panel_i}=THERMAL_output.T{panel_i} +273.15;
        Irr{panel_i}=WEATHER_output.Irr{panel_i};
        if TOOLBOX_input.deviceOptic.exportAbsorptanceAll
            ind = CELL_output.CELL_FRONT.Absmat_ind(2:end-1)-2;
            J_{panel_i}=squeeze(WEATHER_output.J{panel_i}(:,:,ind))*1.6022e-19;
        else
            J_{panel_i}=WEATHER_output.J{panel_i}*1.6022e-19;
        end
    end
end

[TOOLBOX_input] = USR_OPT_ELECTRIC_main(TOOLBOX_input,TYPE); %running the script to obtain user input
Mod_type = TOOLBOX_input.electric.TYPE;
TrackerType = TOOLBOX_input.electric.MPPTrackerType;
minVoltage = TOOLBOX_input.electric.minVoltage;

% Consider Metalization
if(TOOLBOX_input.electric.runMetalization)

    [TOOLBOX_input.electric.shading,TOOLBOX_input.electric.resistance]=...
        metallization(Acell,TOOLBOX_input.electric.numBusbars,...
        TOOLBOX_input.electric.numFingers,TOOLBOX_input.electric.BusbarWdith,TOOLBOX_input.electric.FingerWdith);

end

% Set vaialbes of determining the modules STC power

%Temperatures in Kelvin
T_STC=298.15;
%Conversion from mA/cm2 to A/m2 and into local variable
if TOOLBOX_input.deviceOptic.exportAbsorptanceAll
    J_STC=CELL_output.CELL_FRONT.Jph(1,ind)'*10;
else
    J_STC=CELL_output.CELL_FRONT.Jph(1,:)'*10;
end
T_STC_cells = cell(1,N_submod);
J_STC_cells = cell(1,N_submod);
Irr_STC_cells = cell(1,N_submod);
for SubMod_i = 1:N_submod
    %Set the STC values for each cell
    J_STC_i = J_STC(SUBMOD_IND == SubMod_i);
    T_STC_cells{SubMod_i} = ones(1,numCells(SubMod_i))*T_STC;
    J_STC_cells{SubMod_i} = ones(1,numCells(SubMod_i),length(J_STC_i)).*reshape(J_STC_i,[1,1,length(J_STC_i)]);
    Irr_STC_cells{SubMod_i} = ones(1,numCells(SubMod_i))*1000;
end

for iter = 1:N_iter
    disp('Calculating IV curves for the Modules. This may take a few minutes...');
    % Determine Module IV curve

    if TOOLBOX_input.runPeriodic
        Losses = nan;
        if strcmp(TYPE,'Tan') || strcmp(TYPE,'BIF-Tan')%Tandem modules need a different script
            [V_module,I,Parameters_1, Parameters_2] = TandemModule(TOOLBOX_input,numCells, Acell,T_,J_,J_STC,CELL_output.SUBMOD_IND,Losses);
            [V_module_STC,I_STC,Parameters_STC_1, Parameters_STC_2] = TandemModule(TOOLBOX_input,numCells, Acell,T_STC_cells,J_STC_cells,J_STC,CELL_output.SUBMOD_IND,Losses);
            reconfig_setting = nan;
        elseif strcmp(TYPE,'3Tan') || strcmp(TYPE,'BIF-3Tan')
            [V_module,I,Parameters_1, Parameters_2, Parameters_3] = TripleTandemModule(TOOLBOX_input,numCells, Acell,T_,J_,J_STC,CELL_output.SUBMOD_IND,Losses);
            [V_module_STC,I_STC,Parameters_STC_1, Parameters_STC_2, Parameters_STC_3] = TripleTandemModule(TOOLBOX_input,numCells, Acell,T_STC_cells,J_STC_cells,J_STC,CELL_output.SUBMOD_IND,Losses);
            reconfig_setting = nan;
        else
            [V_module,I,Parameters_1,reconfig_setting] = StandardModule(TOOLBOX_input,numCells, Acell,T_{1},J_{1},J_STC,Irr{1},Losses);
            Parameters_2 = 0;
            [V_module_STC,I_STC,Parameters_STC_1,~] = StandardModule(TOOLBOX_input,numCells, Acell,T_STC_cells{1},J_STC_cells{1},J_STC,Irr_STC_cells{1});
            Parameters_STC_2 =0;
        end
    else %non-periodic simulations
        N_panels = TOOLBOX_input.Scene.N_panels;
        reconfig_setting = cell(N_panels,1);
        V_module = cell(N_panels,1);
        ind_Diode = cell(N_panels,1);
        I = cell(N_panels,1);
        Parameters_1 = cell(N_panels,1);
        Parameters_2 = cell(N_panels,1);
        V_module_STC = cell(N_panels,1);
        I_STC= cell(N_panels,1);
        Parameters_STC_1 = cell(N_panels,1);
        Parameters_STC_2 = cell(N_panels,1);
        if strcmp(TYPE,'Tan') || strcmp(TYPE,'BIF-Tan')%for tandem modules
            for i = 1:N_panels
                temp = T{i}(ind1:ind2,:);
                J1 = J_top{i}(:,:,ind1:ind2);
                J2 = J_bottom{i}(:,:,ind1:ind2);
                J_new_top = converttoELECTRICformat(J1);
                J_new_bottom = converttoELECTRICformat(J2);
                J_tan(:,:,1) = J_new_top{i};
                J_tan(:,:,2) = J_new_bottom{i};
                [V_module{i},I{i}] =...
                    TandemModule(numCells,Acell,temp,J_tan);
                [V_module_STC{i},I_STC{i}] = ...
                    TandemModule(numCells, Acell,T_STC_cells,J_STC_cells);
            end
        else
            for i = 1:N_panels
                Losses = nan;
                [V_module{i},I{i},Parameters_1{i},reconfig_setting{i},ind_Diode{i}] = StandardModule(TOOLBOX_input,numCells,Acell,T_{i},J_{i},J_STC,Irr{i},Losses);
                Parameters_2{i} = 0;
                [V_module_STC{i},I_STC{i},Parameters_STC_1{i},~,~] = StandardModule(TOOLBOX_input,numCells, Acell,T_STC_cells{1},J_STC_cells{1},J_STC,Irr_STC_cells{1},Losses);
                Parameters_STC_2{i} = 0;
            end
        end
    end

    disp('Calculating DC Output from Panels...');
    % Find MPP,Isc,Voc

    if TOOLBOX_input.runPeriodic
        time=length(THERMAL_output.T{1}(:,1));
        [Pmpp, Impp, Vmpp, Isc, Voc]=CreatingModuleDCOutput(V_module,I,time,Mod_type,TrackerType,minVoltage);
        [P_STC, Impp_STC, Vmpp_STC, Isc_STC, Voc_STC]=CreatingModuleDCOutput(V_module_STC, I_STC,1,Mod_type,TrackerType,0);
        ModuleEnergyYield=sum(Pmpp);
    else %non-periodic simulations
        Pmpp = cell(N_panels,1);
        Impp = cell(N_panels,1);
        Vmpp = cell(N_panels,1);
        Isc = cell(N_panels,1);
        Voc = cell(N_panels,1);
        P_STC = cell(N_panels,1);
        Impp_STC = cell(N_panels,1);
        Vmpp_STC = cell(N_panels,1);
        Isc_STC = cell(N_panels,1);
        Voc_STC = cell(N_panels,1);
        ModuleEnergyYield = cell(N_panels,1);

        for i = 1:N_panels
            time = size(T_{i},1);
            [Pmpp{i}, Impp{i},Vmpp{i}, Isc{i}, Voc{i}]=CreatingModuleDCOutput(V_module{i}, I{i},time,Mod_type,TrackerType,minVoltage);
            [P_STC{i}, Impp_STC{i}, Vmpp_STC{i},Isc_STC{i}, Voc_STC{i}]=CreatingModuleDCOutput(V_module_STC{i}, I_STC{i},1,Mod_type,TrackerType,0);
            ModuleEnergyYield{i}=sum(Pmpp{i});
        end
        if TOOLBOX_input.electric.ConnectModules
            %Calculate the operating conditions of every string
            N_strings = length(TOOLBOX_input.electric.Mod_Con);
            for str_i = 1:N_strings
                Mod_sel = TOOLBOX_input.electric.Mod_Con{str_i};
                Imax = max(max(cell2mat(I(Mod_sel))));
                if isfield(TOOLBOX_input.electric,'TaylorParam')
                    I_str=0:(Imax/500):TOOLBOX_input.electric.TaylorParam.I0max;
                else
                    I_str=0:(Imax/500):Imax;
                end
                V_str = zeros(size(V_module{1}));

                %Add the voltages of all modules in the string
                for mod_i = Mod_sel
                    V_additional = max(interp1(I{mod_i},V_module{mod_i}',I_str,'linear','extrap')',0)';
                    if size(V_str,2) == size(V_additional,2)
                        V_str = V_str + V_additional;
                    else
                        V_str = V_str + V_additional';
                    end
                end
                time = size(T_{1},1);
                [Pmpp_str, Impp_str,Vmpp_str, ~, ~]=CreatingModuleDCOutput(V_str, I_str,time,Mod_type,TrackerType,minVoltage);

                %Update the Impp, Pmpp, and Vmpp of the modules in the string
                for mod_i = Mod_sel
                    Impp{mod_i} = Impp_str;
                    for t = 1:length(Vmpp{mod_i}); Vmpp{mod_i}(t) = interp1(I{mod_i},V_module{mod_i}(t,:),Impp_str(t),'linear','extrap')'; end
                    Pmpp{mod_i} = Vmpp{mod_i}.*Impp{mod_i};
                end
            end
        end
    end

    %Calculate dissapating heat
    if TOOLBOX_input.electric.IncludeDissapatingHeat
        if TOOLBOX_input.runPeriodic
            days = find(sum(J_{1},2)>0);
            [HeatGen,ind_Diode] = CalculateGeneratedHeat(I,Impp,Parameters_1,ind_Diode,numCells,days,T_{1},TOOLBOX_input);

            TOOLBOX_input_new = TOOLBOX_input;
            TOOLBOX_input_new.thermal.HeatGen = HeatGen;

            [THERMAL_output_new, ~] = THERMAL_main(TOOLBOX_input_new, CELL_output, MODULE_output, WEATHER_output);


            T_ = cell(1,N_submod);
            for SubMod_i = 1:N_submod
                T_{SubMod_i}=THERMAL_output_new.T{SubMod_i} +273.15;
            end
        else %non-periodic simulation
            HeatGen = cell(N_panels,1);
            for i = 1:N_panels
                days = find(sum(J_{i},2)>0);
                [HeatGen{i},ind_Diode{i}] = CalculateGeneratedHeat(I{i},Impp{i},Parameters_1{i},ind_Diode{i},numCells,days,T_{i},TOOLBOX_input);
            end
            TOOLBOX_input_new = TOOLBOX_input;
            TOOLBOX_input_new.thermal.HeatGen = HeatGen;
            [THERMAL_output_new, ~] = THERMAL_main(TOOLBOX_input_new, CELL_output, MODULE_output, WEATHER_output);
            T_=cellfun(@(x) x+273.15, THERMAL_output_new.T,'UniformOutput',false);

        end

        ELECTRIC_output.HeatGen = HeatGen;
        ELECTRIC_output.ind_Diode = ind_Diode;
        ELECTRIC_output.T_new = T_;

    end
end

ELECTRIC_output.Isc=Isc;
ELECTRIC_output.Impp=Impp;
ELECTRIC_output.Vmpp=Vmpp;
ELECTRIC_output.DCP=Pmpp;
ELECTRIC_output.Voc=Voc;
ELECTRIC_output.P_STC=P_STC;
ELECTRIC_output.Isc_STC=Isc_STC;
ELECTRIC_output.Voc_STC=Voc_STC;
ELECTRIC_output.Impp_STC=Impp_STC;
ELECTRIC_output.Vmpp_STC=Vmpp_STC;
ELECTRIC_output.ModuleEnergyYield_assumingHourlyData...
    = ModuleEnergyYield;
ELECTRIC_output.reconfig_setting = reconfig_setting;
ELECTRIC_output.Parameters_1 = Parameters_1; % added by youri
ELECTRIC_output.Parameters_2 = Parameters_2; % added by youri
ELECTRIC_output.Parameters_STC_1 = Parameters_STC_1; % added by youri
ELECTRIC_output.Parameters_STC_2 = Parameters_STC_2; % added by youri
if strcmp(TYPE,'3Tan') || strcmp(TYPE,'BIF-3Tan')
    ELECTRIC_output.Parameters_3 = Parameters_3; % added by youri
    ELECTRIC_output.Parameters_STC_3 = Parameters_STC_3; % added by youri
end

if ~TOOLBOX_input.script
    TOOLBOX_input.electric.electricplot = true; %true for GUI version
end
if TOOLBOX_input.electric.electricplot
    % plot IV curve
    plot_IVCurve(WEATHER_output.Period, V_module, I,TOOLBOX_input);

    %plot needs to be edited
    % Plots the energy yield and irradiance

    period=WEATHER_output.Period;
    if TOOLBOX_input.runPeriodic
        absorbed_power= WEATHER_output.A*MODULE_output.Amod;
    else
        absorbed_power= WEATHER_output.A;
        Amod = MODULE_output.Panels.Amod;
        absorbed_power = cellfun(@(x) x.*Amod,absorbed_power,...
            'UniformOutput',false);
    end
    plot_EnergyYield(absorbed_power, Pmpp, time,period,TOOLBOX_input);
end

if ~TOOLBOX_input.script
    disp('Electrical calculation finished.DC Output calculated.')
end

end