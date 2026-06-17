%new function to create TOOLBOX_input struct for the script version
%this includes new lines of code for non-periodic simulation + some minor
%edits
%Co-developed by Karthik Ganapathi Subramanian & Malte Vogt
%Source code courtesy of Malte Vogt & Rudi Santbergen, PVMD Group, TU
%Delft, the Netherlands.

% ©All rights reserved.

%% TO USE THE PVMD TOOLBOX SCRIPT version type the command:
%"TB_script('input.mat')" or under which  name (*.mat) you saved the input
%struct created by this file

function TB_input_creator(outputName)

% Add functions to MATLAB path
here = pwd;
addpath(genpath(here));
[src_folder,~,~] = get_folder_structure;
types_folder = fullfile(src_folder, '1_Cell','Types');
addpath(genpath(types_folder));
sim_folder = fullfile(src_folder, '1_Cell','Sim');
addpath(genpath(sim_folder));

disp('Generating Input Struct for TOOLBOX...');
%this needs to be true for the script version to run
TOOLBOX_input.script=true;

% Specify settings file
TOOLBOX_input.settings = 'settings';

% Enable saving results
TOOLBOX_input.save.enable = false;
TOOLBOX_input.save.root_path = append(pwd,'\results'); % Relative/absolute path to root folder
TOOLBOX_input.save.project_name = 'Test';
TOOLBOX_input.save.user_name = 'XYZ';
TOOLBOX_input.save.description = ''; % Optional

% Indicate which parts of the Toolbox should be run
TOOLBOX_input.runDeviceOptic = true;
TOOLBOX_input.runModulePart = true; %running module part
TOOLBOX_input.runWeatherPart = true; %run weather part
TOOLBOX_input.runThermalPart = true; %run thermal part
TOOLBOX_input.runElectricalPart=true; %run electrical part
TOOLBOX_input.runDegradationPart = true; %run degradation part
TOOLBOX_input.runACConversionPart = false; %run conversion part
TOOLBOX_input.runFinancialPart = false; %run financial part
TOOLBOX_input.runLA = false; %run loss analysis part

% Constants
TOOLBOX_input.constants.h = 6.626069E-34; %Plank constant
TOOLBOX_input.constants.c = 299792458; %Speed of light
TOOLBOX_input.constants.T_STC = 298.15; %Temperature at standard test conditions
TOOLBOX_input.constants.q = 1.60217662e-19; %The charge of an electron
TOOLBOX_input.constants.k = 1.380649e-23; %The boltzmann constant
TOOLBOX_input.constants.T_S = 5778; %The temperature of the sun
TOOLBOX_input.constants.AM0 = 1361; %The extraterestial irradiance
TOOLBOX_input.constants.sigma_SB =  5.6697*1e-8; %Stefan Boltzmann constant
TOOLBOX_input.constants.SolidAngleSun = 67.7e-6; % Solid angle from sun to earth(https://en.wikipedia.org/wiki/Solid_angle)

if TOOLBOX_input.runDeviceOptic
    TOOLBOX_input.deviceOptic.runGenPro = true; %set to true for genpro simulations

    % Loading GenPRO file (*.m from types folder), change for new simulations,
    % change if runGenPro is false, to loading stored optical simulatiom (*.mat
    % from Sim folder)
    TOOLBOX_input.deviceOptic.GenProFile = 'ps_tandem_HZB';
    TOOLBOX_input.deviceOptic.exportAbsorptanceAll = false;

    % this is needed to differentiate between tandem, bifacial and normal
    % modules
    if TOOLBOX_input.deviceOptic.runGenPro
        [~, DeviceInfo] = read_device(TOOLBOX_input.deviceOptic.GenProFile);
        celltype = DeviceInfo{2,2};
    else
        load(fullfile(sim_folder,TOOLBOX_input.deviceOptic.GenProFile),'CELL_output');
        celltype = CELL_output.TYPE;
        clear CELL_output;
    end
else
    TOOLBOX_input.deviceOptic = nan;
end

if TOOLBOX_input.runModulePart
    TOOLBOX_input.runPeriodic = true; %SET TO FALSE FOR NON-PERIODIC SIMULATIONS

    if TOOLBOX_input.runPeriodic %inputs for periodic simulations....
        TOOLBOX_input.Scene.loadSimulation = false;
        TOOLBOX_input.Scene.runLUX = false;
        TOOLBOX_input.Scene.runBackwardTracer=true;
        TOOLBOX_input.Scene.tracking = false;
        %load existing file with ray-tracing results
        %SET TO TRUE FOR NEW LUX SIMULATIONS
        TOOLBOX_input.Scene.LuxLoadFile='MonoPERC.mat';  %loading LUX file

        TOOLBOX_input.Scene.module_mounting.ModTilt=0; %module tilt (deg)
        TOOLBOX_input.Scene.module_mounting.ModAzimuth=0; %Module azimuth [deg] (0=S,90=W,180=N,270=E)
        TOOLBOX_input.Scene.module_mounting.ModMountHeight=50; %Module height above ground [cm]

        TOOLBOX_input.Scene.module_mounting.CellRows=12; %Number of cell rows
        TOOLBOX_input.Scene.module_mounting.CellColumns=6; %Number of cell columns
        TOOLBOX_input.Scene.module_mounting.ModThick=0.5; %Module thickness [cm]
        TOOLBOX_input.Scene.module_mounting.CellSpacing=0.3; %Cell spacing [cm]
        TOOLBOX_input.Scene.module_mounting.EdgeSpacing=1; %Edge spacing [cm]
        TOOLBOX_input.Scene.module_mounting.CellLength=15.675; %Cell length [cm]
        TOOLBOX_input.Scene.module_mounting.CellWidth=15.675; %Cell width [cm]

        TOOLBOX_input.Scene.module_mounting.Albedo_eff = 1;
        if TOOLBOX_input.Scene.module_mounting.Albedo_eff
            TOOLBOX_input.Scene.module_mounting.Albedo = 0.2;
        else
            TOOLBOX_input.Scene.module_mounting.Ground_material = 'Concrete';
        end
        if TOOLBOX_input.Scene.tracking
            TOOLBOX_input.Scene.trackingType = 'HSATS'; %Options are 'HSATS','TSATS', 'DATS'
            if strcmp(TOOLBOX_input.Scene.trackingType,'HSATS')
                TOOLBOX_input.Scene.module_mounting.number_rows = 2;
                TOOLBOX_input.Scene.module_mounting.number_columns= 2;
                TOOLBOX_input.Scene.module_mounting.inter_column_distance = 50;
                TOOLBOX_input.Scene.module_mounting.pole_width = 15;
                TOOLBOX_input.Scene.module_mounting.number_modules_side = 2;
                TOOLBOX_input.Scene.module_mounting.azimuth_axis = 0;

            elseif strcmp(TOOLBOX_input.Scene.trackingType,'TSATS')
                TOOLBOX_input.Scene.module_mounting.number_rows = 2;
                TOOLBOX_input.Scene.module_mounting.number_columns= 2;
                TOOLBOX_input.Scene.module_mounting.inter_column_distance = 264;
                TOOLBOX_input.Scene.module_mounting.pole_width = 20;
                TOOLBOX_input.Scene.module_mounting.tilt_axis = 20;
                TOOLBOX_input.Scene.module_mounting.azimuth_axis = 0;

            elseif strcmp(TOOLBOX_input.Scene.trackingType,'DATS')
                TOOLBOX_input.Scene.module_mounting.number_rows = 2;
                TOOLBOX_input.Scene.module_mounting.number_columns= 2;
                TOOLBOX_input.Scene.module_mounting.inter_column_distance = 264;
                TOOLBOX_input.Scene.module_mounting.pole_width = 40;
                TOOLBOX_input.Scene.module_mounting.number_modules_side_x = 2;
                TOOLBOX_input.Scene.module_mounting.number_modules_side_y = 2;

            end
            TOOLBOX_input.Scene.module_mounting.P_tracker = 25; %Watt
            TOOLBOX_input.Scene.module_mounting.eta_tracker = 0.97; %efficiency in [-]
            TOOLBOX_input.Scene.module_mounting.omega_tracker = 0.3; % deg/s

        else
            TOOLBOX_input.Scene.module_mounting.ModSideSpacing=2; %Module side spacing [cm]
            TOOLBOX_input.Scene.module_mounting.ModRowSpacing=800; %Module row spacing [cm]
            if TOOLBOX_input.Scene.runLUX
                TOOLBOX_input.Scene.module_mounting.NRays=500; %25000 is recommended for ray-tracing
                TOOLBOX_input.Scene.ParallelComputation = 1;
                TOOLBOX_input.Scene.MinimumIntensityTermination = 1e-4;
                TOOLBOX_input.Scene.MinimumIntensityMonteCarlo = 0.5;
                TOOLBOX_input.Scene.MaximumGenerationRays = 100;
                TOOLBOX_input.Scene.plot = 0;
                
            else
                TOOLBOX_input.Scene.N_refinement_normal =2;
                TOOLBOX_input.Scene.N_refinement_reduced =1;
            end
        end

        TOOLBOX_input.Scene.plotScene = 0;
        TOOLBOX_input.Scene.module_mounting.avgSensitivity=true;...
            %individual sensitivity for cells
        %SET TO TRUE FOR INDIVIDUAL CELL SENSITIVITY
    else %inputs for non-periodic simulations
        TOOLBOX_input.Scene.N_panels = 4; %number of panels
        TOOLBOX_input.Scene.NonPeriodic_Environment = 'Example_house.mat'; %No environment is selected

        TOOLBOX_input.Scene.loadSimulation = false;
        TOOLBOX_input.Scene.runBackwardTracer=true;
        TOOLBOX_input.Scene.SimulationFile='MonoPERC.mat';
        %load existing file with ray-tracing results
        %SET TO TRUE FOR NEW LUX SIMULATIONS
        TOOLBOX_input.Scene.LuxLoadFile='';

        Tilt = [33.7, 33.7, 33.7, 33.7];
        x_Location = [-250,-150,150,250];
        y_Location = [-200,-200,-200,-200];
        z_Location = [1100/3,1100/3,1100/3,1100/3];
        for Panel_i = 1:TOOLBOX_input.Scene.N_panels
            TOOLBOX_input.Scene.module_mounting(Panel_i).ModTilt=Tilt(Panel_i); %module tilt (deg)
            TOOLBOX_input.Scene.module_mounting(Panel_i).ModAzimuth=0; %Module azimuth [deg] (0=S,90=W,180=N,270=E)

            TOOLBOX_input.Scene.module_mounting(Panel_i).CellRows=12; %Number of cell rows
            TOOLBOX_input.Scene.module_mounting(Panel_i).CellColumns=6; %Number of cell columns
            TOOLBOX_input.Scene.module_mounting(Panel_i).ModThick=0.5; %Module thickness [cm]
            TOOLBOX_input.Scene.module_mounting(Panel_i).CellSpacing=0.3; %Cell spacing [cm]
            TOOLBOX_input.Scene.module_mounting(Panel_i).EdgeSpacing=1; %Edge spacing [cm]
            TOOLBOX_input.Scene.module_mounting(Panel_i).CellLength=15.675; %Cell length [cm]
            TOOLBOX_input.Scene.module_mounting(Panel_i).CellWidth=15.675; %Cell width [cm]
            TOOLBOX_input.Scene.module_mounting(Panel_i).xCoordinate = x_Location(Panel_i); %x Location (bottom left corner) [cm]
            TOOLBOX_input.Scene.module_mounting(Panel_i).yCoordinate = y_Location(Panel_i); %x Location (bottom left corner) [cm]
            TOOLBOX_input.Scene.module_mounting(Panel_i).zCoordinate = z_Location(Panel_i); %x Location (bottom left corner) [cm]
        end

        TOOLBOX_input.Scene.plotScene = 0;
        TOOLBOX_input.Scene.N_refinement_normal =2;
        TOOLBOX_input.Scene.N_refinement_reduced =1;
    end
else
    TOOLBOX_input.Scene = nan;
end

if TOOLBOX_input.runWeatherPart
    % User choice of the plotting of the weather figures
    TOOLBOX_input.irradiation.plot_weather = false;

    %select a file from the locations folder to simulate
    TOOLBOX_input.irradiation.climateFile='Delft-DL-31deg';

    % Select the simulation period
    TOOLBOX_input.irradiation.init_day = 1;
    TOOLBOX_input.irradiation.init_month = 6;
    TOOLBOX_input.irradiation.end_day = 4;
    TOOLBOX_input.irradiation.end_month = 6;

    % Year of simulation
    TOOLBOX_input.irradiation.year_choice = 2021;

    TOOLBOX_input.irradiation.include_horizon_reconstruction = true;

    if TOOLBOX_input.irradiation.include_horizon_reconstruction
        TOOLBOX_input.irradiation.plot_skyline = false;
        TOOLBOX_input.irradiation.skyline_file = 'reconstructed_horizon_DelftMarkt.mat';
    end

    % Choose the spectra that needs to be used:
    %   1 is SMARTS (only consideres clear sky
    %   2 is SBDarts (also considers clouds)
    TOOLBOX_input.irradiation.spectra_choice = 1;
else
    TOOLBOX_input.irradiation = nan;
end

if TOOLBOX_input.runThermalPart
    TOOLBOX_input.thermal.Type = 1; %1 for normal solar cell, 2 for PVT
    TOOLBOX_input.thermal.plot_thermal = false;

    if TOOLBOX_input.thermal.Type == 1
        TOOLBOX_input.thermal.Temperature_model = 'Sandia model';
        if strcmp(TOOLBOX_input.thermal.Temperature_model,'Fluid dynamic model')
            % Efficiency of the PV module [-]
            %(Be aware this is only used for the thermal model for the temperture calculation)
            TOOLBOX_input.thermal.cell_eff = 0.28;
            % Temperature coefficient of the PV module [1/K]
            % cell_eff_temp = cell_eff*(1 + temp_coeff*(module_temp - 300))
            TOOLBOX_input.thermal.temp_coeff = -0.003;
            % Thickness of the glass [m]
            TOOLBOX_input.thermal.glass_thickness = 0.0032;
            TOOLBOX_input.thermal.anenometer_height = 10; %[m]
            TOOLBOX_input.thermal.glass_conduct = 1; %[W/(m K)]
            TOOLBOX_input.thermal.glass_emissivity = 0.89; %[-]
            
            %The thermal response is estimated at 7 minutes, so the
            %following value is advised: exp(-T_step [minutes]/ 7 [minutes])
            TOOLBOX_input.thermal.T_tau = exp(-60/7);

        elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Duffie-Beckman model')
            % Efficiency of the PV module [-]
            %(Be aware this is only used for the thermal model for the temperture calculation)
            TOOLBOX_input.thermal.cell_eff = 0.28;

            % NOCT temperature
            TOOLBOX_input.thermal.T_NOCT = 60;

        elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Faiman model')
            % Efficiency of the PV module [-]
            %(Be aware this is only used for the thermal model for the temperture calculation)
            TOOLBOX_input.thermal.cell_eff = 0.28;

            % constant heat transfer component
            TOOLBOX_input.thermal.U0 = 25;

            % convective heat transfer component
            TOOLBOX_input.thermal.U1 = 6.84;

            % alpha
            TOOLBOX_input.thermal.alpha = 0.9;

        elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Sandia model')
            % Parameters for the Sandia model
            TOOLBOX_input.thermal.a = -3.47;
            TOOLBOX_input.thermal.b = -0.0594;

        elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Incropera model')
            %The number of layers
            TOOLBOX_input.thermal.Nlayers = 5;

            %Enable rear convection (1 = yes, 0 = no);
            TOOLBOX_input.thermal.RearConvection = 1;

            % Fixed temperature at the rear (NaN means no fixed temperature)
            TOOLBOX_input.thermal.RearTemperature = nan;

            % Efficiency of the PV module [-]
            %(Be aware this is only used for the thermal model for the temperture calculation)
            TOOLBOX_input.thermal.Efficiency = 0.20;


            %Set up the structure for the Incropera model
            TOOLBOX_input.thermal.layers  = [PvLayer('name','Enc top','thickness',1/1000,'k',0.24,'rho',1700,'cp',250,'dx',0.25,'nIntNodeY',6,'emissivity',0.89)
                PvLayer('name','EVA top','thickness',0.46*2/1000,'k',0.32,'rho',960,'cp',2090,'dx',0.25,'nIntNodeY',1)
                PvLayer('name','Cell','thickness',0.2/1000,'k',149,'rho',2330,'cp',838,'dx',0.25,'nIntNodeY',3)
                PvLayer('name','EVA bottom','thickness',0.46/1000,'k',0.32,'rho',960,'cp',2090,'dx',0.25,'nIntNodeY',1)
                PvLayer('name','Enc bot','thickness',0.4/1000,'k',0.56,'rho',1370,'cp',1760,'dx',0.25,'nIntNodeY',1,'emissivity',0.89)];

            % Assign to which layer all materials belong
            TOOLBOX_input.thermal.assignment_layers = [0,1,1,1,2,2,2,2,3,3,3,3,3,0];
        end
    elseif TOOLBOX_input.thermal.Type == 2
        TOOLBOX_input.thermal.runpvt = 'Solar_PV_thermal_I_tank.m';
    end
else
    TOOLBOX_input.thermal = nan;
end

if TOOLBOX_input.runElectricalPart
    TOOLBOX_input.electric.InterpolateIV =false; %Use the interpolation method

    if TOOLBOX_input.electric.InterpolateIV
        %The parameters correspond to:
        %[Rs0 [Ohm], Rp0 [Ohm], Is1,0 [A], n1 [-], Is2,0 [A], n2 [-], Eg [eV],
        %T-Iph [-], TR-p1 [-], TR-s1 [-], TX-Is1 [-], TX-Is2 [-]]
        TOOLBOX_input.electric.IVparameters = [0.0022, 192.8402, 3.9961e-12, 1, 5.7427e-07, 2, 1.14, 2.9e-3, 0, 0, 3, 3];

        %The IV conditions correspond to:
        %[Min irradiance [W/m^2],step irradiance [W/m^2], max irradiance [W/m^2], min temperature [C], step temperature [C], max temperature [C]]
        TOOLBOX_input.electric.IVconditions = [0, 50, 1200, -10, 5,100];

        TOOLBOX_input.electric.PVmodtype = 'sp';
    else %MODEL CELL VALUES
        if strcmp(celltype,'Tan')|| strcmp(celltype,'BIF-Tan')
            TOOLBOX_input.electric.IVtypeTop = 'PerovskiteV2024';
            TOOLBOX_input.electric.IVtypeBot = 'SiliconV2024';
            TOOLBOX_input.electric.LC_eff = 0;
            TOOLBOX_input.electric.degtype = [2,1];
            TOOLBOX_input.electric.k_needed = {0,0};
            TOOLBOX_input.electric.degScenario = {1,nan};

            TOOLBOX_input.electric.Terminals= 2; %can be 2, 3, or 4 terminal tandem
            if TOOLBOX_input.electric.Terminals == 3
                TOOLBOX_input.electric.VM_ratio_m = 3;
                TOOLBOX_input.electric.VM_ratio_n = 2;
            end


        elseif strcmp(celltype,'3Tan') || strcmp(celltype,'BIF-3Tan')
            TOOLBOX_input.electric.IVtypeTop = 'PerovskiteV2024';
            TOOLBOX_input.electric.IVtypeMid = 'PerovskiteV2024';
            TOOLBOX_input.electric.IVtypeBot = 'SiliconV2024';
            TOOLBOX_input.electric.LC_eff = [0,0,0];
            TOOLBOX_input.electric.degtype = [2,2,1];
            TOOLBOX_input.electric.k_needed = {0,0,0};
            TOOLBOX_input.electric.degScenario = {1,1,nan};
            TOOLBOX_input.electric.tripleType = 'ABC'; %Options are 'ABC', 'A-BC', or 'AB-C'
        else
            %Applies to  single junction modules
            TOOLBOX_input.electric.IVtype = 'SiliconV2024';
            TOOLBOX_input.electric.degtype = 1;
            TOOLBOX_input.electric.k_needed = {0};
            TOOLBOX_input.electric.degScenario = {nan};
        end
    end

    %Set degradation inputs to defaul as they are not used in the
    %electrical part yet
    TOOLBOX_input.electric.deg_sil.k_mois = 0;
    TOOLBOX_input.electric.deg_sil.factor_LID = 1;
    TOOLBOX_input.electric.deg_sil.filename = 'TrendMoistureReference';

    TOOLBOX_input.electric.runMetalization=false; %metalization not-considered

    if TOOLBOX_input.electric.runMetalization
        TOOLBOX_input.electric.numBusbars=5; %Number of busbars (#/cell)
        TOOLBOX_input.electric.numFingers=130; %Number of fingers (#/cell)
        TOOLBOX_input.electric.BusbarWdith=500e-6; %Busbar thickness (m)
        TOOLBOX_input.electric.FingerWdith=30e-6;%Finger thickness(m)
    else
        TOOLBOX_input.electric.shading=2; %specify Loss due to shading by metalization[%]
        TOOLBOX_input.electric.resistance=0.0058; %Resistance of cell metal grid [Ohm]
    end

    TOOLBOX_input.electric.numBypassDiodes=3; %bypass diodes can be 1,3 or 6
    TOOLBOX_input.electric.forwardVoltageDiode = 0;
    TOOLBOX_input.electric.electricplot = false; %no plots for script version
    TOOLBOX_input.electric.ModuleType = 'Series'; %for reconfigurable modules, change this value to 'REC'
    TOOLBOX_input.electric.MPPTrackerType = 'Global';
    TOOLBOX_input.electric.minVoltage = 0;
    TOOLBOX_input.electric.ConnectModules = 1;
    TOOLBOX_input.electric.Mod_Con = {[1,2,3,4]}; %The strings should be indicated with arrays (i.e. {[string 1],[string 2]}
    TOOLBOX_input.electric.IncludeDissapatingHeat = 0;
    TOOLBOX_input.electric.N_iter = 2;
else
    TOOLBOX_input.electric = nan;
end

if TOOLBOX_input.runDegradationPart == true
    TOOLBOX_input.Degradation.loadCOMSOL = true;
    TOOLBOX_input.Degradation.N_years = 20; %Amount of simulated years
    TOOLBOX_input.Degradation.COMSOL_file = 'RMC_Delft';
    TOOLBOX_input.Degradation.MoistureDeg_file = 'TrendMoistureReference';
    TOOLBOX_input.Degradation.Discoloration_file = 'TrendDiscolorationReference';
    if strcmp(celltype,'Tan')|| strcmp(celltype,'BIF-Tan')
        TOOLBOX_input.Degradation.DegScenario = 1;
        TOOLBOX_input.Degradation.FactorDegPerov = 1;
    end
else %no degradation calculation
    TOOLBOX_input.Degradation = nan;
end

if TOOLBOX_input.runACConversionPart
    TOOLBOX_input.Conversion.plot_conversion = false;

    % Type of inverter. Options: 'Central', 'String', 'Micro', 'Power opt.'
    TOOLBOX_input.Conversion.ConversionType = 'String';
    TOOLBOX_input.Conversion.Model='ABB: PVI 3.0 OUTD-S-US-Z-M-A (208 V) 208V [CEC 2014]';
    if strcmp(TOOLBOX_input.Conversion.Model,'Specify inverter')
        TOOLBOX_input.Conversion.inv_constants(1) = nan;
        TOOLBOX_input.Conversion.inv_constants(2) = nan;
        TOOLBOX_input.Conversion.inv_constants(3) = nan;
        TOOLBOX_input.Conversion.inv_constants(4) = 12000;
        TOOLBOX_input.Conversion.inv_constants(5) = 12269.9387;
        TOOLBOX_input.Conversion.inv_constants(6) = 600;
        TOOLBOX_input.Conversion.inv_constants(7) = 10.5269;
        TOOLBOX_input.Conversion.inv_constants(8) = -3.85E-07;
        TOOLBOX_input.Conversion.inv_constants(9) = -1.62E-05;
        TOOLBOX_input.Conversion.inv_constants(10) = -0.0014;
        TOOLBOX_input.Conversion.inv_constants(11) = -0.0073;
        TOOLBOX_input.Conversion.inv_constants(12) = 0;
        TOOLBOX_input.Conversion.inv_constants(13) = 1100;
    else
        TOOLBOX_input.Conversion.inv_constants = obtain_inverter_constants(TOOLBOX_input.Conversion.Model);
    end

    if strcmp(TOOLBOX_input.Conversion.ConversionType,'Power opt.')
        TOOLBOX_input.Conversion.PowerOpt.Model = 'P-300';
        TOOLBOX_input.Conversion.PowerOpt.AddCentralInverter = true;
        if TOOLBOX_input.Conversion.PowerOpt.AddCentralInverter
            TOOLBOX_input.Conversion.Model = 'ABB: PVI 3.0 OUTD-S-US-Z-M-A (208 V) 208V [CEC 2014]';
            % Fixed voltage only applicable for periodic simulations
            TOOLBOX_input.Conversion.PowerOpt.FixedVoltage = false;
            % Determine the number of modules in parallel and in series
            TOOLBOX_input.Conversion.Parallel_Modules = 2;
            TOOLBOX_input.Conversion.Series_Modules = 2;
        end
    elseif strcmp(TOOLBOX_input.Conversion.ConversionType,'Micro')
        % The variables below are only used to know the number of PV modules in
        % the system, being Parallel_Modules*Series_Modules
        TOOLBOX_input.Conversion.Parallel_Modules=1;
        TOOLBOX_input.Conversion.Series_Modules=1;
    elseif strcmp(TOOLBOX_input.Conversion.ConversionType,'String')
        TOOLBOX_input.Conversion.Series_Modules=6;
        % The variable below stands for the number of strings in the system
        TOOLBOX_input.Conversion.Parallel_Modules=1;
    else
        TOOLBOX_input.Conversion.Parallel_Modules=2;
        TOOLBOX_input.Conversion.Series_Modules=4;
    end

    %Cable losses: specify if you want to perform a detailed calculation or use
    %a fixed percentage
    TOOLBOX_input.Conversion.CableLossesDetailedCalculation = false;
    if TOOLBOX_input.Conversion.CableLossesDetailedCalculation
        % Electrical resistivity in ohm*mm^2/m 0.0168 for copper, 0.02820 for
        % aluminium
        TOOLBOX_input.Conversion.InverterCableResistivity = 0.0168;
        % Length of the cable in m
        TOOLBOX_input.Conversion.InverterCableLength = 20;
        % Cross section of the cable in mm2
        TOOLBOX_input.Conversion.InverterCableCrossSection = 4;

        if contains(['Central','Power opt.'],TOOLBOX_input.Conversion.ConversionType)
            % Only for the central inverter and for the power optimizer when
            % there is a central inverter, one can specify the cable losses of
            % the string from the modules to the junction box
            TOOLBOX_input.Conversion.StringCableResistivity = 0.0168;
            TOOLBOX_input.Conversion.StringCableLength = 6;
            TOOLBOX_input.Conversion.StringCableCrossSection = 4;
        end
    else
        %Specify the fixed percentage for cable losses
        TOOLBOX_input.Conversion.CableLoss = 0.5;
    end
else
    TOOLBOX_input.Conversion = nan;
end

if TOOLBOX_input.runFinancialPart == true
    if TOOLBOX_input.irradiation.starting_month == 1 && ...
            TOOLBOX_input.irradiation.starting_day == 1 && ...
            TOOLBOX_input.irradiation.ending_month == 12 && ....
            TOOLBOX_input.irradiation.ending_day == 31  %for yearly simulations
        %include code here for LCOE calculations
        TOOLBOX_input.FinancialPart.DiscountRate = 5;
        TOOLBOX_input.FinancialPart.Modulecost = 395.47;
        TOOLBOX_input.FinancialPart.Invertercost = 151.45;
        TOOLBOX_input.FinancialPart.StructuralBOScost = 84.14;
        TOOLBOX_input.FinancialPart.ElectricalBOScost = 193.53;
        TOOLBOX_input.FinancialPart.OMcost = 16.8;
        TOOLBOX_input.FinancialPart.Inverterlifetime = 50;
        TOOLBOX_input.FinancialPart.Installyear = 13;
        TOOLBOX_input.FinancialPart.scenario = 1; %default scenario for LCOE
    else
        TOOLBOX_input.runFinancialPart = false;
        TOOLBOX_input.FinancialPart = nan; %no financial data included...
    end
else %no LCOE calculation
    TOOLBOX_input.FinancialPart = nan; %no financial data included...
end

if TOOLBOX_input.runLA
    TOOLBOX_input.LossAnalysis.Run_Operating = 1;
    TOOLBOX_input.LossAnalysis.TrackingType = 0;
    TOOLBOX_input.LossAnalysis.plotFigures = 0;
    TOOLBOX_input.LossAnalysis.B_STC = 0; % set 1 for considering rear side irradiance
    if strcmp(celltype,'Tan') || strcmp(celltype,'BIF-Tan')
        TOOLBOX_input.LossAnalysis.E_g1 = 1.68;
        TOOLBOX_input.LossAnalysis.E_g2 = 1.12;
    else
        TOOLBOX_input.LossAnalysis.E_g = 1.12;
    end
else
    TOOLBOX_input.LossAnalysis = nan;
end

if exist('outputName', 'var')==0
    outputName='input';
end

save(outputName, 'TOOLBOX_input'); %save variable
disp('Input Struct Defined Succesfully!');
disp('Struct saved as:');
disp([outputName,'.mat']);
end


function inv_constants = obtain_inverter_constants(inverterModel)
%GET_INVERTER_CONSTANTS Return SNL constants of the inverter selected
% Read the inverters characteristics on SNL database and select those that
% are compatible: if max array voltage does not exceed max inverter input
% voltage. Compare that selected database with the inverter chosen by the
% user
%
% Parameters
% ----------
% Vdc_max : double
%   Maximum DC voltage in the time period
% inverter : struct
%   User inputs related to the inverter characteristics
%
% Returns
% -------
% inv_constants : double
%   SNL constants of the chosen inverter if compatible

% Import data from SNL database
[constants, models] = xlsread('InverterConstants.xlsx');
models = models(4:end,1);


model_idx = find(strcmp(models,inverterModel),1);
if isempty(find(strcmp(models,inverterModel),1))
    error("conversion:InvalidInverterModel", "The selected inverter model %s is not in the list",...
        inverter.model);
else
    inv_constants = constants(model_idx,:);
end

end