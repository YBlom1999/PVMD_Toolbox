function Losses_STC = Analysis_STC(TOOLBOX_input,CELL_output,MODULE_output,ELECTRIC_output,CONVERSION_output)
%Analysis_Operating Calculates the losses at STC.
%
% This function calculates the losses that are present in the PV system for
% STC
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Simulation results of the CELL module
% MODULE_output : struct
%   Simulation results of the MODULE module
% ELECTRIC_output : struct
%   Simulation results of the ELECTRIC module
% CONVERSION_output : struct
%   Simulation results of the CONVERSION module
%
% Returns
% -------
% Losses_STC : struct
%   Simulation results of the loss analysis at STC
%
% Developed by Y. Blom

B_STC = 0; % 1 for B-STC simulation, keep 0 as default.
[src_folder,~,~] = get_folder_structure;
load(fullfile(src_folder,'1_CELL','GenPro4','spectrum.mat'),'spec'); %The AM1.5 spectrum is loaded

%% a structure of all constants are defined
CONSTANTS.h = 6.62607004e-34;
CONSTANTS.q = 1.60217662e-19;
CONSTANTS.c = 299792458;
CONSTANTS.k = 1.380649e-23;
CONSTANTS.T_S = 5778;

T_cell = 298.15;

%Properties of the module
A_mod = MODULE_output.Amod;
N_cells = MODULE_output.N;
if B_STC == 0
    Angle_emit = 2*pi;
elseif B_STC == 1
    Angle_emit = 4*pi;
end
Angle_abs = 67.7e-6; %https://en.wikipedia.org/wiki/Solid_angle
V_mpp = ELECTRIC_output.Vmpp_STC;
I_mpp = ELECTRIC_output.Impp_STC;
type = CELL_output.TYPE;

if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F') %If it is a single junction
    Parameters = ELECTRIC_output.Parameters_STC_1(1,:);
    Terminals = 2;
    E_g = TOOLBOX_input.LossAnalysis.E_g;
    R = CELL_output.CELL_FRONT.RAT(:,1,1)+CELL_output.CELL_FRONT.RAT(:,1,end);
    A_diff = 1 - sum(CELL_output.CELL_FRONT.RAT(:,1,:),3);
    A = CELL_output.CELL_FRONT.RAT(:,1,2);
elseif strcmp(type,'Tan') || strcmp(type,'BIF-Tan') %If it is a tandem juntion
    Parameters1 = ELECTRIC_output.Parameters_STC_1(1,:);
    Parameters2 = ELECTRIC_output.Parameters_STC_2(1,:);
    Terminals = TOOLBOX_input.electric.Terminals;
    E_g1 = TOOLBOX_input.LossAnalysis.E_g1;
    E_g2 = TOOLBOX_input.LossAnalysis.E_g2;
    R = CELL_output.CELL_FRONT.RAT(:,1,1)+CELL_output.CELL_FRONT.RAT(:,1,end);
    A_diff = 1 - sum(CELL_output.CELL_FRONT.RAT(:,1,:),3);
    A1 = CELL_output.CELL_FRONT.RAT(:,1,2);
    A2 = CELL_output.CELL_FRONT.RAT(:,1,3);
end



%% Data

%The spectrum is loaded and transferred into the spectral irradiance and
%the spectral photon flux
h = CONSTANTS.h;
c = CONSTANTS.c;
q = CONSTANTS.q;
wav = spec.data(:,1)*1e-6;
Irr_spec = spec.data(:,2)*1e10*h*c./(q*wav);
photon_spec = spec.data(:,2)*1e10/q;
for i = 1 : length(wav)-1
    Irr_spec(i) = Irr_spec(i)*1e-9/(wav(i+1)-wav(i));
    photon_spec(i) = photon_spec(i)*1e-9/(wav(i+1)-wav(i));
end
if B_STC == 1
    Irr_spec = 1.2*Irr_spec;
    photon_spec = 1.2*photon_spec;
end


%% Power in
%A numberical integral is performed to calculate the total irradiance power
P_in = A_mod*trapz(wav,Irr_spec);

% The optimal voltage is calculated for the given irradiance level
if strcmp(type,'SHJ') || strcmp(type,'BIF')
    V_opt1 = Vopt_calculator(wav,photon_spec',E_g,mean(T_cell),Angle_emit);
elseif strcmp(type,'Tan') || strcmp(type,'BIF-Tan')
    V_opt1 = Vopt_calculator(wav,photon_spec',E_g1,mean(T_cell),Angle_emit);
    V_opt2 = Vopt_calculator(wav,photon_spec',E_g2,mean(T_cell),Angle_emit);
end

%% Fundamental losses
if strcmp(type,'SHJ') || strcmp(type,'BIF')|| strcmp(type,'T-F')
    Fund_Losses_input.wav = wav;
    Fund_Losses_input.Irr_spec = Irr_spec;
    Fund_Losses_input.photon_spec = photon_spec;
    Fund_Losses_input.Angle_abs = Angle_abs;
    Fund_Losses_input.Angle_emit_Vopt = Angle_emit;
    Fund_Losses_input.Angle_emit_emission = Angle_emit;
    Fund_Losses_input.A = A;
    Fund_Losses_input.E_g = E_g;
    Fund_Losses_input.I_mpp = I_mpp;
    Fund_Losses_input.V_mpp = V_mpp;
    Fund_Losses_input.T_cell = T_cell;
    Fund_Losses_input.V_opt1 = V_opt1;
    [P_term,P_below,P_emission,P_carnot,P_angle,P_gain] = FundamentalLossesSingle(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);
elseif strcmp(type,'Tan') || strcmp(type,'BIF-Tan')
    Fund_Losses_input.wav = wav;
    Fund_Losses_input.Irr_spec = Irr_spec;
    Fund_Losses_input.photon_spec = photon_spec;
    Fund_Losses_input.Angle_abs = Angle_abs;
    Fund_Losses_input.Angle_emit_Vopt = Angle_emit;
    Fund_Losses_input.Angle_emit_emission = Angle_emit;
    Fund_Losses_input.A1 = A1;
    Fund_Losses_input.A2 = A2;
    Fund_Losses_input.E_g1 = E_g1;
    Fund_Losses_input.E_g2 = E_g2;
    Fund_Losses_input.Parameters1 = Parameters1;
    Fund_Losses_input.I_mpp = I_mpp;
    Fund_Losses_input.V_mpp = V_mpp;
    Fund_Losses_input.T_cell = T_cell;
    Fund_Losses_input.V_opt1 = V_opt1;
    Fund_Losses_input.V_opt2 = V_opt2;
    [P_term,P_below,P_emission,P_carnot,P_angle,P_gain] = FundamentalLossesTandem(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);
end
P_fund = P_term+P_below+P_emission+P_carnot+P_angle-P_gain;

%% Optical losses
if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F')
    Opt_Losses_input.P_in = P_in;
    Opt_Losses_input.P_fund = P_fund;
    Opt_Losses_input.wav = wav;
    Opt_Losses_input.photon_spec = photon_spec;
    Opt_Losses_input.R = R;
    Opt_Losses_input.A = A;
    Opt_Losses_input.A_diff = A_diff;
    Opt_Losses_input.E_g = E_g;
    Opt_Losses_input.V_opt1 = V_opt1;
    Opt_Losses_input.Angle_emit = Angle_emit;
    Opt_Losses_input.V_mpp = V_mpp;
    Opt_Losses_input.I_mpp = I_mpp;
    Opt_Losses_input.T_cell = T_cell;
    [P_cell,P_metal,P_ref,P_diffA,I_abs] = OpticalLossesSingle(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);
elseif strcmp(type,'Tan') || strcmp(type,'BIF-Tan')
    Opt_Losses_input.P_in = P_in;
    Opt_Losses_input.P_fund = P_fund;
    Opt_Losses_input.wav = wav;
    Opt_Losses_input.photon_spec = photon_spec;
    Opt_Losses_input.R = R;
    Opt_Losses_input.A1 = A1;
    Opt_Losses_input.A2 = A2;
    Opt_Losses_input.A_diff = A_diff;
    Opt_Losses_input.E_g1 = E_g1;
    Opt_Losses_input.E_g2 = E_g2;
    Opt_Losses_input.V_opt1 = V_opt1;
    Opt_Losses_input.V_opt2 = V_opt2;
    Opt_Losses_input.Angle_emit = Angle_emit;
    Opt_Losses_input.Parameters1 = Parameters1;
    Opt_Losses_input.V_mpp = V_mpp;
    Opt_Losses_input.I_mpp = I_mpp;
    Opt_Losses_input.T_cell = T_cell;
    [P_cell,P_metal,P_ref,P_diffA,I_abs1,I_abs2] = OpticalLossesTandem(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);
end

if strcmp(type,'SHJ') || strcmp(type,'BIF')
    IV_curvename = TOOLBOX_input.electric.IVtype;
    [P_diffA, I_abs] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs,T_cell,V_opt1',IV_curvename);
elseif strcmp(type,'Tan')|| strcmp(type,'BIF-Tan')    
    IV_curvename_top = TOOLBOX_input.electric.IVtypeTop;
    IV_curvename_bot = TOOLBOX_input.electric.IVtypeBot;
    [P_diffA, I_abs1] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs1,T_cell,V_opt1',IV_curvename_top);
    [P_diffA, I_abs2] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs2,T_cell,V_opt2',IV_curvename_bot);
end


%% Electrical losses
if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F')
    Elec_Losses_input.I_abs = I_abs;
    Elec_Losses_input.V_opt1 = V_opt1;
    Elec_Losses_input.Parameters = ones(N_cells,1,5).*reshape(Parameters,1,1,5);

    Elec_Losses_input.T_cell = T_cell;
    [P_series,P_shunt,P_NRRI,P_NRRV] = ElectricLossesSingle_avg(Elec_Losses_input,MODULE_output);
elseif strcmp(type,'Tan') || strcmp(type,'BIF-Tan')
    Elec_Losses_input.I_abs1 = I_abs1;
    Elec_Losses_input.I_abs2 = I_abs2;
    Elec_Losses_input.V_opt1 = V_opt1;
    Elec_Losses_input.V_opt2 = V_opt2;
    Elec_Losses_input.Parameters1 = ones(N_cells,1,5).*reshape(Parameters1,1,1,5);
    Elec_Losses_input.Parameters2 = ones(N_cells,1,5).*reshape(Parameters2,1,1,5);
    Elec_Losses_input.T_cell = T_cell;
    [P_series,P_shunt,P_NRRI,P_NRRV,PowerRatio] = ElectricLossesTandem_avg(Elec_Losses_input,TOOLBOX_input,MODULE_output,CONSTANTS);
end

%% System losses
if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F')
    Sys_Losses_input.I_mpp = I_mpp;
    Sys_Losses_input.V_mpp = V_mpp;
    Sys_Losses_input.I_abs = I_abs;
    Sys_Losses_input.Parameters = ones(N_cells,1,5).*reshape(Parameters,1,1,5);
    Sys_Losses_input.T_cell = T_cell;
    if (isfield(TOOLBOX_input, 'runACConversionPart')==1) %To check whether the AC simulation is performed
        if TOOLBOX_input.runACConversionPart == 1
            Sys_Losses_input.Pdc = ELECTRIC_output.P_STC;
            Sys_Losses_input.Pac = CONVERSION_output.Pac_STC;
        end
    end
    [P_con,P_mismatch,P_cable,P_inv] = SystemLossesSingle(Sys_Losses_input,TOOLBOX_input,ELECTRIC_output,MODULE_output);
elseif strcmp(type,'Tan') || strcmp(type,'BIF-Tan')
    if Terminals == 2 || Terminals ==3
        Sys_Losses_input.I_mpp = I_mpp;
        Sys_Losses_input.V_mpp = V_mpp;
        Sys_Losses_input.I_abs1 = I_abs1;
        Sys_Losses_input.I_abs2 = I_abs2;
        Sys_Losses_input.Parameters1 = ones(N_cells,1,5).*reshape(Parameters1,1,1,5);
        Sys_Losses_input.Parameters2 = ones(N_cells,1,5).*reshape(Parameters2,1,1,5);
        Sys_Losses_input.T_cell = T_cell;

        if (isfield(TOOLBOX_input, 'runACConversionPart')==1 && TOOLBOX_input.runACConversionPart == 1) %To check whether the AC simulation is performed
            Sys_Losses_input.Pdc = ELECTRIC_output.P_STC;
            Sys_Losses_input.Pac = CONVERSION_output.Pac_STC;
        end

        [P_con,P_mismatch,P_cable,P_inv] = SystemLossesTandem2T(Sys_Losses_input,TOOLBOX_input,ELECTRIC_output,MODULE_output);
    end
end

%% Total power
disp('DC losses')
P_out = ELECTRIC_output.P_STC;
Power = [P_term;P_below;P_angle;P_NRRV;P_diffA;P_carnot;P_cell;P_emission;P_NRRI;P_ref;P_metal;P_con;P_mismatch;P_shunt;P_series;-P_gain;P_out];
P_total = sum(Power);
Power = [Power; P_total];
Components = ["Thermalization";"Below bandgap";"Angle mismatch";"Recombination V";"Parasitic absorption";"Carnot losses";"Cell spacing";"Emission losses";"Recombination I";"Reflection/Transmission";"Metal shading";"Cell interconnection";"Mismatch losses";"Shunt resistance";"Series resistance";"Gain";"Power";"Total"];
Losses_STC.Power_DC = Power;
Percentage = 100*Power/P_in;
disp(table(Components,Power,Percentage))
Losses_STC.Components_DC = Components;
Losses_STC.Percentage_DC = Percentage;

if (isfield(TOOLBOX_input, 'runACConversionPart')==1)
    if TOOLBOX_input.runACConversionPart == 1
        disp('AC losses')
        Pac_STC = CONVERSION_output.Pac_STC;
        panels = TOOLBOX_input.Conversion.Parallel_Modules*TOOLBOX_input.Conversion.Series_Modules;
        Power = panels*[P_term;P_below;P_angle;P_NRRV;P_diffA;P_carnot;P_cell;P_emission;P_NRRI;P_inv/panels;P_ref;P_metal;P_cable/panels;P_con;P_mismatch;P_shunt;P_series;-P_gain;Pac_STC/panels];
        P_total = sum(Power);
        Power = [Power; P_total];
        Components = ["Thermalization";"Below bandgap";"Angle mismatch";"Recombination V";"Parasitic absorption";"Carnot losses";"Cell spacing";"Emission losses";"Recombination I";"Inverter losses";"Reflection/Transmission";"Metal shading";"Cable losses";"Cell interconnection";"Mismatch losses";"Shunt resistance";"Series resistance";"Gain";"Power";"Total"];
        Losses_STC.Power_AC = Power;
        Percentage = 100*Power/(P_in*panels);
        disp(table(Components,Power,Percentage))
        Losses_STC.Components_AC = Components;
        Losses_STC.Percentage_AC = Percentage;
    end
end

%% Plot results
if TOOLBOX_input.LossAnalysis.plotFigures == 1
    if (isfield(TOOLBOX_input, 'runACConversionPart')==1 && TOOLBOX_input.runACConversionPart == 1)
        if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F') %AC simulation for single junction
            Plot_Components = ["AC power";"Inverter losses";"Cable losses";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Categories = ["Fundamental losses:", "Optical losses:", "Electrical losses:","System losses:","Output Power"];
            Plot_Power = 100*[Pac_STC;P_inv;P_cable;P_mismatch*panels;P_con*panels;P_NRRI*panels;P_NRRV*panels;P_shunt*panels;P_series*panels;P_diffA*panels;P_ref*panels;P_metal*panels;P_cell*panels;-P_gain*panels;P_angle*panels;P_carnot*panels;P_emission*panels;P_below*panels;P_term*panels]/(P_in*panels);
            Plot_Power_Categories = 100*[P_term*panels+P_below*panels+P_angle*panels+P_carnot*panels+P_emission*panels-P_gain*panels;P_cell*panels+P_metal*panels+P_ref*panels+P_diffA*panels;P_NRRV*panels+P_NRRI*panels+P_series*panels+P_shunt*panels;P_mismatch*panels+P_con*panels+P_cable+P_inv;Pac_STC]/(P_in*panels);
            Plot_type = 1;
        elseif strcmp(type,'Tan')  || strcmp(type,'BIF-Tan')  %AC simulation for tandem cell
            Plot_Components = ["AC power bottom cell";"AC power top cell";"Inverter losses";"Cable losses";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Categories = ["Fundamental losses:", "Optical losses:", "Electrical losses:","System losses:","Output Power"];
            Plot_Power = 100*[Pac_STC*(1-PowerRatio);Pac_STC*PowerRatio;P_inv;P_cable;P_mismatch*panels;P_con*panels;P_NRRI*panels;P_NRRV*panels;P_shunt*panels;P_series*panels;P_diffA*panels;P_ref*panels;P_metal*panels;P_cell*panels;-P_gain*panels;P_angle*panels;P_carnot*panels;P_emission*panels;P_below*panels;P_term*panels]/(P_in*panels);
            Plot_Power_Categories = 100*[P_term*panels+P_below*panels+P_angle*panels+P_carnot*panels+P_emission*panels-P_gain*panels;P_cell*panels+P_metal*panels+P_ref*panels+P_diffA*panels;P_NRRV*panels+P_NRRI*panels+P_series*panels+P_shunt*panels;P_mismatch*panels+P_con*panels+P_cable+P_inv;Pac_STC]/(P_in*panels);
            Plot_type = 2;
        end
    else
        if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F') %DC simulation for single junction
            Plot_Components = ["DC power";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Categories = ["Fundamental losses:", "Optical losses:", "Electrical losses:","System losses:","Output Power"];
            Plot_Power = 100*[P_out;P_mismatch;P_con;P_NRRI;P_NRRV;P_shunt;P_series;P_diffA;P_ref;P_metal;P_cell;-P_gain;P_angle;P_carnot;P_emission;P_below;P_term]/P_in;
            Plot_Power_Categories = 100*[P_term+P_below+P_angle+P_carnot+P_emission-P_gain;P_cell+P_metal+P_ref+P_diffA;P_NRRV+P_NRRI+P_series+P_shunt;P_mismatch+P_con;P_out]/P_in;
            Plot_type = 3;
        elseif strcmp(type,'Tan')  || strcmp(type,'BIF-Tan')%DC simulation for tandem cell
            Plot_Components = ["DC power bottom cell";"DC power top cell";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Categories = ["Fundamental losses:", "Optical losses:", "Electrical losses:","System losses:","Output Power"];
            Plot_Power = 100*[P_out*(1-PowerRatio);Pdc_STC*PowerRatio;P_mismatch;P_con;P_NRRI;P_NRRV;P_shunt;P_series;P_diffA;P_ref;P_metal;P_cell;-P_gain;P_angle;P_carnot;P_emission;P_below;P_term]/P_in;
            Plot_Power_Categories = 100*[P_term+P_below+P_angle+P_carnot+P_