function [MODULE_output, TOOLBOX_input] = MODULE_nonperiodic(TOOLBOX_input, CELL_output)
%WEATHER_periodic Calculates the sensitivity map for non-periodic simulations
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
% Developed by Youri Blom


if TOOLBOX_input.script==false
    TOOLBOX_input.runPeriodic = false; %setting periodic simulation as false
    Q ={'Select modelling option'};    %ask user
    choice = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Load saved simulation','Perform backward ray tracing'}); %user choice
    if isempty(choice), return, end
    if choice ==1
        [src_folder,~,~] = get_folder_structure;
        list = dir(fullfile(src_folder,'2_MODULE','Non-Periodic','Simulations',CELL_output.TYPE,'*.mat'));              %list the file names there
        Q ='Choose the Simulation:';    %ask user
        A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
        if isempty(A), return, end
        
        %copy choices to workspace
        TOOLBOX_input.Scene.loadSimulation=true;
        TOOLBOX_input.Scene.runLUX=false;
        TOOLBOX_input.Scene.runBackwardTracer=false;
        TOOLBOX_input.Scene.SimulationFile=list(A).name;
        
    elseif choice == 2
        %set run lux
        TOOLBOX_input.Scene.loadSimulation=false;
        TOOLBOX_input.Scene.runBackwardTracer=true;
        
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        cd(fullfile(pvmd_folder, '2_MODULE','Non-Periodic','Environments')); %go to 'Environment' folder (where the environments are stored)
        list = dir('*.mat');              %list the file names there
        Q ='Choose the Environment:';    %ask user
        A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
        if isempty(A), return, end
        
        N_panels = input('Specify number of modules: ');
        %copy choices to workspace
        TOOLBOX_input.Scene.NonPeriodic_Environment = list(A).name;
        cd(current_path);
        TOOLBOX_input.Scene.N_panels = N_panels;
    end
end
    
if TOOLBOX_input.Scene.loadSimulation
    [src_folder,~,~] = get_folder_structure;
    load(fullfile(src_folder,'2_MODULE','Non-Periodic','Simulations',CELL_output.TYPE,TOOLBOX_input.Scene.SimulationFile),'MODULE_output')
    % TOOLBOX_input_mount = load(fullfile(src_folder,'2_MODULE','Non-Periodic','Simulations',CELL_output.TYPE,TOOLBOX_input.Scene.SimulationFile),'TOOLBOX_input');
    % if exist('TOOLBOX_input_mount','var')
    %     TOOLBOX_input.Scene.module_mounting = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting;
    %     TOOLBOX_input.Scene.N_panels = TOOLBOX_input_mount.TOOLBOX_input.Scene.N_panels;
    % 
    % end
    TOOLBOX_input.Scene.N_panels = length(MODULE_output.SM_f);
else
    %Load the environment
    [src_folder,~,~] = get_folder_structure;
    load(fullfile(src_folder,'2_MODULE','Non-Periodic','Environments',TOOLBOX_input.Scene.NonPeriodic_Environment),'V','F','materials','opa','rgb','Albedo','Scattering')
    opa_env = opa;
    rgb_env = rgb;
    N_env = size(F,1);
    Albedo_env = Albedo;
    Scattering_env = Scattering;
    V = 100*V;
    Panels.cellCorners = [];
    Panels.normalSolarCell = [];
    Panels.Ncells = [];
    Panels.materials = [];
    Panels.T = [];
    Panels.Acell = [];
    Panels.Amod = [];
    Panels = repelem(Panels,1,TOOLBOX_input.Scene.N_panels);
    for Panel_i = 1:TOOLBOX_input.Scene.N_panels
        [V, F,materials_i,BF,cellCorners_i,normalSolarCell_i,Ncells_i,Acell_i,Amod_i,TOOLBOX_input,T_i] = AddPanels_NonPeriodic(TOOLBOX_input,CELL_output,V,F,materials,Panel_i);
        if BF
            Panels(Panel_i).cellCorners = cellCorners_i{1};
            Panels(Panel_i).cellCorners_rear = cellCorners_i{2};
        else
            Panels(Panel_i).cellCorners = cellCorners_i;
        end
        Panels(Panel_i).normalSolarCell = normalSolarCell_i;
        Panels(Panel_i).Ncells = Ncells_i;
        Panels(Panel_i).materials = materials_i;
        Panels(Panel_i).T = T_i;
        Panels(Panel_i).Acell = Acell_i;
        Panels(Panel_i).Amod = Amod_i;
    end
    
    % plot_Environment(F,V,N_env,rgb_env,opa_env,TOOLBOX_input,Panels);
    
    
    %===create an array for the albedo and scattering of every vertex
    wav = CELL_output.CELL_FRONT.wav;
    Albedo = zeros(length(F),length(wav));
    Scattering = zeros(1,length(F));
    Albedo(1:N_env,:) = Albedo_env.*ones(N_env,length(wav));
    Scattering(1:N_env) = Scattering_env;
    for Panel_i = 1:TOOLBOX_input.Scene.N_panels
        T = Panels(Panel_i).T;
        for t = 1:length(T)
            if ~isempty(T(t).Facet)
                index_start = 2*T(t).Facet(1)-1;
                index_end = 2*(T(t).Facet(end)-T(t).Facet(1))+index_start+1;
                index_rgb = (index_start:index_end)+N_env+2*T(end).Facet(end)*(Panel_i-1);
                if isstruct(T(t).RT)
                    Albedo(index_rgb,:) = mean(T(t).RT.RAT(:,:,1),2)'.*ones(length(index_rgb),length(wav));
                else
                    Albedo(index_rgb,:) = T(t).RT(1);
                end
                if isfield(T(t),'Scat') && ~isempty(T(t).Scat)
                    Scattering(index_rgb) = 1;
                end
            end
        end
    end
    
    
    %Calculate SM map
    N_panels = TOOLBOX_input.Scene.N_panels;
    SM_f = cell(N_panels,1);
    if BF; SM_r = cell(N_panels,1); end
    
    for Panel_i = 1:TOOLBOX_input.Scene.N_panels
        Ncells = Panels(Panel_i).Ncells;
        cellCorners = Panels(Panel_i).cellCorners;
        normalSolarCell = Panels(Panel_i).normalSolarCell;
        T = Panels(Panel_i).T;
        [SM_f{Panel_i},Vs,Fs,azimuth,zenith,As] = BackwardTracer(V,F,Ncells,cellCorners,normalSolarCell,Albedo,Scattering,T(1).RT);
        % flatplot3(Vs,Fs,mean(mean(SM_f{Panel_i}(:,:,1,:),4),2),Panel_i+1);
        if BF
            cellCorners_rear = Panels(Panel_i).cellCorners_rear;
            [SM_r{Panel_i},~,~,~,~,~] = BackwardTracer(V,F,Ncells,cellCorners_rear,-normalSolarCell,Albedo,Scattering,T(2).RT);
            % plot_i = TOOLBOX_input.Scene.N_panels+Panel_i;
            % flatplot3(Vs,Fs,mean(mean(SM_r{Panel_i}(:,:,1,:),4),2),plot_i);
        end
    end
    
    MODULE_output.skydome.AZA = [azimuth,zenith,As];    %skydome zenith, azimuth, area for every triangle center
    MODULE_output.skydome.Vs = Vs;                      %skydome vertices (for plotting SKYmap)
    MODULE_output.skydome.Fs = Fs;                      %skydome facets (for plotting SKYmap)
    MODULE_output.SM_f = SM_f;
    MODULE_output.wav = CELL_output.CELL_FRONT.wav;     %pass on wavelength information
    MODULE_output.Panels = Panels;
    MODULE_output.ModTilt = TOOLBOX_input.Scene.module_mounting.ModTilt;
    if BF, MODULE_output.SM_r = SM_r; end
end
end

function plot_Environment(F,V,N_env,rgb_env,opa_env,TOOLBOX_input,Panels)
%plot_Environment plots the environment of the simulation
%
% This function makes a figure of the environment and the PV modules
%
% Parameters
% ----------
% F : double
%   The faces of the environment (including the modules)
% V : double
%   The vertices in the environment (including the modules)
% N_env : double
%   The number of faces belonging to the environment (excluding the modules)
% rgb_env : double
%   The color of faces belonging to the environment (excluding the modules)
% opa_env : double
%   The opacity of faces belonging to the environment (excluding the modules)
% TOOLBOX_input : struct
%   Simulation parameters
% Panels : struct
%   Properties of all the PV panels
%
% Returns
% -------
%
% Developed by Youri Blom
nrf = size(F,1);

rgb = zeros(nrf,3);                      %initialize
opa = zeros(nrf,1);                      %initialize
rgb(1:N_env,:) = rgb_env;
opa(1:N_env) = opa_env;

for Panel_i = 1:TOOLBOX_input.Scene.N_panels
    T = Panels(Panel_i).T;
    for t = 1:length(T)
        if ~isempty(T(t).Facet)
            index_start = 2*T(t).Facet(1)-1;
            index_end = 2*(T(t).Facet(end)-T(t).Facet(1))+index_start+1;
            index_rgb = (index_start:index_end)+N_env+2*T(end).Facet(end)*(Panel_i-1);
            if isstruct(T(t).RT)
                rgb(index_rgb,3) = 1;      %give cell blue color
                opa(index_rgb) = 1-mean(mean(T(t).RT.RAT(:,:,end)));
            else
                rgb(index_rgb,:) = T(t).RT(1)*ones(2*numel(T(t).Facet),3);   %others gray-scale
                opa(index_rgb) = 1-T(t).RT(2);
            end
        end
    end
end
p = patch('Faces',F,'Vertices',V);       %plot facets
%set rgb and opacity
set(p,'FaceVertexCData',rgb,'CDataMapping','scaled',...
    'FaceColor','flat','FaceVertexAlphaData',opa,...
    'AlphaDataMapping','none','FaceAlpha','flat')

axis equal off tight                    %xyz aspect ratio equal
view(30,30)                            %standard view from side
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-1000,1000]);
ylim([-1000,1000]);
zlim([0,1000]);

end