function [TOOLBOX_input] = USR_OPT_ELECTRIC_main(TOOLBOX_input,TYPE)
% Select Modelling options
if TOOLBOX_input.script==false
    Q ={'Select modelling option'};    %ask user
    Model = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Datasheet','Model Cells','Interpolation of IV curves'},'ListSize',[150,100]); %user choice
    if isempty(Model), return, end    %if user presses 'cancel'
    if Model ==1 %Datasheet
        TOOLBOX_input.electric.runDatasheet=true;
        TOOLBOX_input.electric.InterpolateIV =false;
        
        %Determine datasheet values
        if strcmp(TYPE,'Tan') || strcmp(TYPE,'BIF-Tan')
            Q = {'Voc [V]','Isc [A]','Vmpp [V]','Impp [A]','n [-]','Kp [% / ?C]','Kv [% / ?C]','Ki [% / ?C]'};
            default = {'93.6','4.55','83.3','4.32','1.3','-0.35','-0.30','0.03'};
            TOOLBOX_input.electric.datasheetValuesTop = str2double(inputdlg(Q,'Top Module',1,default));
            
            Q = {'Voc [V]','Isc [A]','Vmpp [V]','Impp [A]','n [-]','Kp [% / ?C]','Kv [% / ?C]','Ki [% / ?C]'};
            default = {'47.52','4.55','38.5','4.32','1.3','-0.35','-0.30','0.03'};
            TOOLBOX_input.electric.datasheetValuesBot = str2double(inputdlg(Q,'Bottom Module',1,default));
        elseif strcmp(TYPE,'3Tan') || strcmp(TYPE,'BIF-3Tan')
            Q = {'Voc [V]','Isc [A]','Vmpp [V]','Impp [A]','n [-]','Kp [% / ?C]','Kv [% / ?C]','Ki [% / ?C]'};
            default = {'93.6','4.55','83.3','4.32','1.3','-0.35','-0.30','0.03'};
            TOOLBOX_input.electric.datasheetValuesTop = str2double(inputdlg(Q,'Top Module',1,default));
            
            Q = {'Voc [V]','Isc [A]','Vmpp [V]','Impp [A]','n [-]','Kp [% / ?C]','Kv [% / ?C]','Ki [% / ?C]'};
            default = {'93.6','4.55','83.3','4.32','1.3','-0.35','-0.30','0.03'};
            TOOLBOX_input.electric.datasheetValuesMiddle = str2double(inputdlg(Q,'Top Module',1,default));
            
            Q = {'Voc [V]','Isc [A]','Vmpp [V]','Impp [A]','n [-]','Kp [% / ?C]','Kv [% / ?C]','Ki [% / ?C]'};
            default = {'47.52','4.55','38.5','4.32','1.3','-0.35','-0.30','0.03'};
            TOOLBOX_input.electric.datasheetValuesBot = str2double(inputdlg(Q,'Bottom Module',1,default));

        else
            Q = {'Voc [V]','Isc [A]','Vmpp [V]','Impp [A]','n [-]','Kp [% / ?C]','Kv [% / ?C]','Ki [% / ?C]'};
            default = {'40','8','35','7.5','1.3','-0.35','-0.30','0.03'};
            TOOLBOX_input.electric.datasheetValues = str2double(inputdlg(Q,'Values',1,default));
        end
        
        
    elseif Model ==2 %I-V fittings
        TOOLBOX_input.electric.runDatasheet=false;
        TOOLBOX_input.electric.InterpolateIV =false;
        
        %Select IV curve
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        cd(fullfile(pvmd_folder, '5_ELECTRIC','Model','Data')); %go to 'Data' folder
        list = dir('*.mat');              %list the file names there
        cd(current_path);
        if strcmp(TYPE,'Tan') || strcmp(TYPE,'BIF-Tan')
            Q ='Select the cell IV curve type for the top cells:';
            A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
            if isempty(A), return, end
            
            TOOLBOX_input.electric.IVtypeTop = list(A).name;
            
            
            Q ='Select the cell IV curve type for the bot cells:';
            A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
            if isempty(A), return, end
            
            TOOLBOX_input.electric.IVtypeBot = list(A).name;
            TOOLBOX_input.electric.factors_top = [1,1,1,1,1];
            TOOLBOX_input.electric.factors_bot = [1,1,1,1,1];
        elseif strcmp(TYPE,'3Tan') || strcmp(TYPE,'BIF-3Tan')
            Q ='Select the cell IV curve type for the top cells:';
            A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
            if isempty(A), return, end
            
            TOOLBOX_input.electric.IVtypeTop = list(A).name;

            Q ='Select the cell IV curve type for the middle cells:';
            A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
            if isempty(A), return, end
            
            TOOLBOX_input.electric.IVtypeMid = list(A).name;
            
            
            Q ='Select the cell IV curve type for the bot cells:';
            A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
            if isempty(A), return, end
            
            TOOLBOX_input.electric.IVtypeBot = list(A).name;
        else
            Q ='Select the cell IV curve type for the cells:';
            A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
            if isempty(A), return, end
            
            TOOLBOX_input.electric.IVtype = list(A).name;
            TOOLBOX_input.electric.factors = [1,1,1,1,1];
        end
    elseif Model == 3 % Interpolation of IV curves
        TOOLBOX_input.electric.runDatasheet=false;
        TOOLBOX_input.electric.InterpolateIV =true;
        % Ask for the cell parameters
        Q = {'Rs0 [Ohm]','Rp0 [Ohm]','Is1,0 [A]','n1 [-]','Is2,0 [A]','n2 [-]','Eg [eV]','T-Iph [-]', 'TR-p1 [-]','TR-s1 [-]','TX-Is1 [-]','TX-Is2 [-]'};
        default = {'0.0022','192.8402','3.9961e-12','1','5.7427e-07','2','1.14','2.9e-3','0','0','3','3'};
        TOOLBOX_input.electric.IVparameters = str2double(inputdlg(Q,'Values',1,default));
        
        % Ask for the different conditions
        Q = {'minimum irradiance [W/m^2]','step irradiance [W/m^2]', 'maximum irradiance [W/m^2]','minimum temperature [C]', 'step temperature [C]','maximum temperature [C]'};
        default = {'0','50','1200','-10','5','100'};
        TOOLBOX_input.electric.IVconditions = str2double(inputdlg(Q,'Values',1,default));
        
        Q ={'Select module type'};    %ask user
        PVmodtype = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'series-parallel','butterfly'},'ListSize',[150,100]); %user choice
        switch PVmodtype
            case 1
                TOOLBOX_input.electric.PVmodtype = 'sp';
            case 2
                TOOLBOX_input.electric.PVmodtype = 'but';
        end
    end
    
    
    
    
    ch = input('Are the Modules Re-configurable? [Y/N]: ', 's');
    if strcmp(ch,'Y')
        numCells = TOOLBOX_input.Scene.module_mounting(1).CellRows*TOOLBOX_input.Scene.module_mounting(1).CellColumns;
        if numCells ~= 96 %if 6 blocks not possible
            warning('Re-configurable Modules only possible for a Module with 96 cells');
            disp('Module Type Overridden');
            TOOLBOX_input.electric.TYPE = 'Non-REC';
        elseif numCells == 96
            nc_r = TOOLBOX_input.Scene.module_mounting.CellRows;
            nc_c = TOOLBOX_input.Scene.module_mounting.CellColumns;
            if nc_r == 8 && nc_c ==12 || nc_r == 12 && nc_c ==8
                TOOLBOX_input.electric.TYPE = 'REC';
                prompt1 = {'Enter Switch Resistance (m-ohm)',...
                    'Loss due to shading by metalization[%]',...
                    'Resistance of a single cell in a block [Ohm]'};
                dlgtitle = 'Resistance data for Re-configurable Modules';
                dims = [1 40];
                definput = {'2','3','0.00263'};
                R_rec = inputdlg(prompt1,dlgtitle,dims,definput);
                Rsw = str2double(R_rec{1})*10^-3;
                TOOLBOX_input.electric.shading= str2double(R_rec{2});
                TOOLBOX_input.electric.resistance = str2double(R_rec{3});
                TOOLBOX_input.electric.runMetalization=false;
                algo_options = {'SCCS','OSD+IR','OSD+SCCS','OSD+IRwoAP','OSD+SCCSwoAP',...
                    'CC4AP+IR','CC4AP+SCCS','NCC+IR','NCC+SCCS'};
                algo = listdlg('PromptString','Select a reconfiguration algorithm',...
                    'SelectionMode','single','ListString',algo_options);
                TOOLBOX_input.electric.reconfig.algo_nm = algo_options{algo};
                TOOLBOX_input.electric.reconfig.algo = algo;
                TOOLBOX_input.electric.reconfig.Rsw = Rsw;
                TOOLBOX_input.electric.reconfig.nc_r = nc_r;
            else
                warning('Re-configurable Module possible only for Cell Arrangements 8x12 & 12x8!');
                disp('Module Type Overridden');
                TOOLBOX_input.electric.TYPE = 'Non-REC';
            end
        end
    else
        TOOLBOX_input.electric.TYPE = 'Non-REC';
        
        Q ={'How do you whish to consider metalization losses?'};
        Metalization = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Ignore metalization losses','Give shading and resistance values','Detailed model for metalization (Developed for square shaped c-Si cells, unvalidated for other cells)'},'ListSize',[600,100]); %user choice
        if isempty(Model), return, end    %if user presses 'cancel'
        if Metalization ==1
            TOOLBOX_input.electric.shading=0;
            TOOLBOX_input.electric.resistance=0;
            TOOLBOX_input.electric.runMetalization=false;
        elseif Metalization ==2
            if strcmp(TOOLBOX_input.electric.TYPE,'Non-REC')
                Q = {'Loss due to shading by metalization[%]','Resistance of cell metal grid [Ohm]'};
                default = {'3','0.0051'};
                A = str2double(inputdlg(Q,'Metallization',1,default));
                TOOLBOX_input.electric.shading=A(1);
                TOOLBOX_input.electric.resistance=A(2);
            else
                Q = input('Enter the Resistance of cell metal grid [Ohm]: ');
                TOOLBOX_input.electric.shading = Q;
            end
            TOOLBOX_input.electric.runMetalization=false;
        elseif Metalization ==3
            Q = {'Number of busbars (#/cell)','Number of fingers (#/cell)','Busbar thickness (m)','Finger thickness(m)'};
            default = {'5','130','500e-6','30e-6'};
            A = str2double(inputdlg(Q,'Metallization',1,default));
            TOOLBOX_input.electric.numBusbars=A(1);
            TOOLBOX_input.electric.numFingers=A(2);
            TOOLBOX_input.electric.BusbarWdith=A(3);
            TOOLBOX_input.electric.FingerWdith=A(4);
            TOOLBOX_input.electric.runMetalization=true;
        end
        
        %number of bypass diodes
        Q ='Select the number of bypass diodes:';
        List={'0','3','6'};
        Bypassdiodes=[1 3 6];
        TOOLBOX_input.electric.numBypassDiodes = ...
            Bypassdiodes(listdlg('PromptString',Q,...
            'SelectionMode','single','ListString',List));
    end
    if strcmp(TYPE,'Tan') || strcmp(TYPE,'BIF-Tan')
        %Tandem Configuration
        Q ='Select the tandem configuration:'; %ask user
        List={'2 Terminal','3 Terminal','4 Terminal'}; %changed by Youri
        Terminals=[2 3 4]; %changed by Youri
        TOOLBOX_input.electric.Terminals=Terminals...
            (listdlg('PromptString',Q,...
            'SelectionMode','single','ListString',...
            List,'ListSize',[200,100])); %number of terminals
        
        if TOOLBOX_input.electric.Terminals == 3
            TOOLBOX_input.electric.configuration = questdlg('Which configuration?','Configuration choice','Series','Reverse','Reverse');
            TOOLBOX_input.electric.VM_ratio_m = input('Which value for m (top cells)? ');
            TOOLBOX_input.electric.VM_ratio_n = input('Which value for n (bottom cells)? ');
        end
    end
    if strcmp(TYPE,'3Tan') || strcmp(TYPE,'BIF-3Tan')
        %Tandem Type
        Q ='Select the tandem type:'; %ask user
        List={'A (top and middle+bottom)','B (top+middle and bottom','C (top+middle+bottom'}; %changed by Youri
        Types=['A','B','C']; %changed by Youri
        TOOLBOX_input.electric.Type=Types(listdlg('PromptString',Q,...
            'SelectionMode','single','ListString',...
            List,'ListSize',[200,100])); %number of terminals

        % Ask for luminescence coupling efficiencies
        Q = {'LC efficiency top-> mid [%]','LC efficiency top-> bot [%]','LC efficiency mid-> bot [%]'};
        default = {'0','0','0'};
        TOOLBOX_input.electric.LC_eff = str2double(inputdlg(Q,'Metallization',1,default));
    end
end
end