function [MODULE_output, TOOLBOX_input] = MODULE_periodic(TOOLBOX_input, CELL_output)
%WEATHER_periodic Calculates the sensitivity map for periodic simulations
%
% This function calculates the sensitivity map of the PV system
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Output of the CELL block
%
% Returns
% -------
% MODULE_output : struct
%   Output of this module block
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by Karthik Ganapathi Subramanian


%load input values
TOOLBOX_input.runPeriodic = true; %setting periodic simulation as true
TOOLBOX_input.Scene.NonPeriodic_Environment = '';

if TOOLBOX_input.script==false
    Q ={'Select modelling option'};    %ask user
    choice = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Load saved simulation','Perform forward ray tracing','Perform backward ray tracing'}); %user choice
    if isempty(choice), return, end
    if choice ==1
        [src_folder,~,~] = get_folder_structure;
        list = dir(fullfile(src_folder,'2_MODULE','Periodic','Simulations',CELL_output.TYPE,'*.mat'));              %list the file names there
        Q ='Choose the cell technology:';    %ask user
        A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
        if isempty(A), return, end
        
        %copy choices to workspace
        TOOLBOX_input.Scene.loadSimulation=true;
        TOOLBOX_input.Scene.runLUX=false;
        TOOLBOX_input.Scene.runBackwardTracer=false;
        TOOLBOX_input.Scene.LuxLoadFile=list(A).name;
    elseif choice == 2
        %set run lux
        TOOLBOX_input.Scene.loadSimulation=false;
        TOOLBOX_input.Scene.runLUX=true;
        TOOLBOX_input.Scene.runBackwardTracer=false;
    elseif choice == 3
        %set run backward tracer
        TOOLBOX_input.Scene.loadSimulation=false;
        TOOLBOX_input.Scene.runLUX=false;
        TOOLBOX_input.Scene.runBackwardTracer=true;
    end
end


if TOOLBOX_input.Scene.loadSimulation
    [src_folder,~,~] = get_folder_structure;
    load(fullfile(src_folder,'2_MODULE','Periodic','Simulations',CELL_output.TYPE,TOOLBOX_input.Scene.LuxLoadFile),'MODULE_output')
    TOOLBOX_input_mount = load(fullfile(src_folder,'2_MODULE','Periodic','Simulations',CELL_output.TYPE,TOOLBOX_input.Scene.LuxLoadFile),'TOOLBOX_input');
    if exist('TOOLBOX_input_mount','var')
        TOOLBOX_input.Scene.module_mounting.CellRows = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.CellRows;
        TOOLBOX_input.Scene.module_mounting.CellColumns = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.CellColumns;
        TOOLBOX_input.Scene.module_mounting.ModThick = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.ModThick;
        TOOLBOX_input.Scene.module_mounting.CellSpacing = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.CellSpacing;
        TOOLBOX_input.Scene.module_mounting.EdgeSpacing = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.EdgeSpacing;
        TOOLBOX_input.Scene.module_mounting.ModTilt = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.ModTilt;
        TOOLBOX_input.Scene.module_mounting.ModAzimuth = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.ModAzimuth;
        TOOLBOX_input.Scene.module_mounting.ModMountHeight = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.ModMountHeight;
        TOOLBOX_input.Scene.module_mounting.ModSideSpacing = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.ModSideSpacing;
        TOOLBOX_input.Scene.module_mounting.ModRowSpacing = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.ModRowSpacing;
        TOOLBOX_input.Scene.module_mounting.CellLength = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.CellLength;
        TOOLBOX_input.Scene.module_mounting.CellWidth = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.CellWidth;
        TOOLBOX_input.Scene.module_mounting.Albedo = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting.Albedo;
        
    end
    if TOOLBOX_input.script==false
        %===ask user input about module geometry===
        Q = {'Module tilt [deg]',...
            'Module azimuth [deg] (0=S,90=W,180=N,270=E)',...
            'Height to ground [cm]'};
        
        default = {'27','0','50'};
        
        GEO = inputdlg(Q,'MODULE MOUNTING INFORMATION NEEDED FOR THERMAL MODEL',1,default);
        if isempty(GEO), return, end      %if user presses 'cancel'
        
        
        %Module tilt [deg]
        TOOLBOX_input.Scene.module_mounting.ModTilt=str2double(GEO{1});
        
        %Module azimuth [deg] (0=S,90=W,180=N,270=E)
        TOOLBOX_input.Scene.module_mounting.ModAzimuth=str2double(GEO{2});
        
        %Module height above ground [cm]
        TOOLBOX_input.Scene.module_mounting.ModMountHeight=str2double(GEO{3});
    end
elseif TOOLBOX_input.Scene.runLUX
    [MODULE_output, TOOLBOX_input] = runLUX_main(TOOLBOX_input, CELL_output);
elseif TOOLBOX_input.Scene.runBackwardTracer
    [MODULE_output, TOOLBOX_input] = runBackwardTracer_main(TOOLBOX_input, CELL_output);
end

end
