function [CELL_output,TOOLBOX_input] = CELL_main(TOOLBOX_input)
% This function is called when the CELL button is pressed.
%
% The user selects the cell type and the output is the reflectance,
% absorptance in every layer and transmittance as a function of wavelength
% and angle of incidence, both for front and rear-side illumination.
%
% Returns
% -------
% struct
%   Cell_output with fields:
%
%       - RAT_f: front absorptance (dim1: wavelengths, dim2: angles, dim3: layers)
%       - RAT_r: rear absorptance
%       - wav: list of wavelenghts
%       - aoi: list of angles of inc
%       - lay: list of layer names
%


%===Select cell type and get corresponding structure (Lay and Int)===
if TOOLBOX_input.script==false
    Q ={'Select modelling option:'};    %ask user
    choice = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Load saved simulation','Run GenPro4'},'ListSize',[150,80]); %user choice
    if isempty(choice), return, end
    if choice==1
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        cd(fullfile(pvmd_folder, '1_CELL','Sim')); %go to 'Sim' folder (where simulations are stored)
        list = dir('*.mat');              %list the file names there
        Q ='Choose the cell technology:';    %ask user
        A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
        if isempty(A), return, end
        
        %copy choices to workspace
        TOOLBOX_input.deviceOptic.runGenPro=false;
        TOOLBOX_input.deviceOptic.GenProFile=list(A).name;
        cd(current_path);
    else
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        cd(fullfile(pvmd_folder, '1_CELL','Types'));    %go to 'Types' folder (where cell types are stored)
        list = dir('*.m');                              %list the file names there
        Q ='Define solar cell type';                    %ask user
        A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name}); %user choice
        if isempty(A), return, end                      %if user presses 'cancel'
        Q = 'Do you want to export the absorptance in all layers (needed for the incropera temperature model)?';
        exportAbsorptanceAll = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'no','yes'},'ListSize',[500,100])-1; %user choice
        if isempty(exportAbsorptanceAll), return, end                      %if user presses 'cancel'
        %[Lay,Int,rback,TYPE] = feval(list(A).name(1:end-2));   %evaluate the selected file (get Lay and Int)
        
        %copy choices to workspace
        TOOLBOX_input.deviceOptic.runGenPro=true;
        TOOLBOX_input.deviceOptic.GenProFile=list(A).name;
        TOOLBOX_input.deviceOptic.exportAbsorptanceAll=exportAbsorptanceAll;
        cd(current_path)                        %return to parent folder
    end
end

if TOOLBOX_input.deviceOptic.runGenPro==false
    current_path = pwd;
    [pvmd_folder,~,~] = get_folder_structure;
    cd(fullfile(pvmd_folder, '1_CELL','Sim')); %go to 'Sim' folder (where simulations are stored)
    load(TOOLBOX_input.deviceOptic.GenProFile,'CELL_output');
    if isfield(CELL_output.CELL_FRONT,'Absmat_ind')
        TOOLBOX_input.deviceOptic.exportAbsorptanceAll = 1;
    else
        TOOLBOX_input.deviceOptic.exportAbsorptanceAll = 0;
    end
    cd(current_path)                        %return to parent folder
else
    %get info from seperately defined GenProFile
    current_path = pwd;
    [pvmd_folder,~,~] = get_folder_structure;
    cd(fullfile(pvmd_folder, '1_CELL','Types'));    %go to 'Types' folder (where cell types are stored)
    [Lay,Int,rback,TYPE,absmat,Submod_ind] = feval(TOOLBOX_input.deviceOptic.GenProFile(1:end-2));
    cd(current_path)
    
    %===Use GP4 to get R,A,T as a function of wavelength and angle===
    disp('CELL calculation started. This can take a few minutes...')
    cd(fullfile(pvmd_folder, '1_CELL','GenPro4'));
    
    %---set angles of incidence---
    S = gp4_settings;                          %load GenPro4 settings
%     S.nai = 15;                                %set nr of angular intervals
    
    if isempty(rback)                   %if rback is empty, bificial module is used
        S.sd = 0;                               %disables an acceleration trick that doesn't work for rear side illumination
    end
    
    [Lay,Int,out] = GENPRO4(Lay,Int,S);    %run Genpro
    title(['Front: ',num2str(0.5*90/S.nai),'^\circ'])
    drawnow
    
    %%%%Absorbermaterials are defined here
    %absmat = {'perovskite','c-Si','perov_2','c-Si-2015', 'perovskite-1.63eV-Man18'};            %recognised absorber materials (TODO: user should be able to modify)
    ix = 1;                                    %include first column (R)
    for c = 1:length(absmat)
        ix = [ix,find(strcmp(out.leg,absmat{c}))]; %#ok<AGROW> %include absorber layers column (A)
    end
    ix = [ix,size(out.abp,1)];                 %include final column (T)
    
    %if the absorptance for all materials is needed, ix is changed
    if TOOLBOX_input.deviceOptic.exportAbsorptanceAll 
        AbsorberMaterials = ix;
        ix = 1:size(out.abp,1);
    end
    
    Jph_STC(1,:) = out.cur(ix(2:end-1));         %[mA/cm2]
    
    RAT_f = zeros(size(out.abp,2),S.nai,length(ix));
    RAT_f(:,1,:) = out.abp(ix,:).'; %save output (first (R), absorber layers (A) and final (T))
    
    
    %---INCIDENCE FROM THE FRONT---
    for c = 2:S.nai       %for every other angular interval on front side
        S.iai = c;        %set incident angular interval
        [~,~,out] = GENPRO4(Lay,Int,S); %run
        title(['Front: ',num2str((c-0.5)*90/S.nai),'^\circ'])
        drawnow
        RAT_f(:,c,:) = out.abp(ix,:).';           %put absorptance output in array
        Jph_STC(c,:) = out.cur(ix(2:end-1));         %[mA/cm2]
    end
    %---put output in struct and copy to workspace---
    CELL_FRONT.RAT = RAT_f;
    CELL_FRONT.wav = S.wav';
    CELL_FRONT.aoi = ((1:S.nai)'-0.5)*90/S.nai;
    CELL_FRONT.lay = out.leg(ix)';
    CELL_FRONT.Jph = Jph_STC;
    %If all absorptance is needed, the index of the absorber materials is
    %needed.
    if TOOLBOX_input.deviceOptic.exportAbsorptanceAll
        CELL_FRONT.Absmat_ind = AbsorberMaterials;
    end
    
    
    %---INCIDENCE FROM THE REAR---
    if isempty(rback)               %if rback is empty, bificial module is used
        RAT_r = zeros(size(RAT_f)); %initialize same matrix as for front side
        Jph_STC_r = zeros(size(Jph_STC));
        for c = 1:S.nai             %for every angular interval on rear side
            S.iai =  S.nai * (4 * length(Int) - 2) + c;      %set incident angular interval (counting back from final one)
            [~,~,out] = GENPRO4(Lay,Int,S); %run
            title(['Rear: ',num2str((c-0.5)*90/S.nai),'^\circ'])
            drawnow
            RAT_r(:,c,:) = out.abp(ix,:).';   %put output data in array
            Jph_STC_r(c,:) = out.cur(ix(2:end-1));         %[mA/cm2]
            %CELL_REAR.Siai_test(c)=S.iai;
        end
        %---put output in struct and copy to workspace---
        CELL_REAR.RAT = RAT_r;
        CELL_REAR.wav = S.wav';
        CELL_REAR.aoi = ((1:S.nai)'-0.5)*90/S.nai;
        CELL_REAR.lay = out.leg(ix)';
        CELL_REAR.Jph = Jph_STC_r;
        if TOOLBOX_input.deviceOptic.exportAbsorptanceAll
            CELL_REAR.Absmat_ind = AbsorberMaterials;
        end
    else                            %otherwise use fixed reflectance value rback
        CELL_REAR = [rback,0];
    end
    
    CELL_output.CELL_FRONT = CELL_FRONT;    %put output in single struct
    CELL_output.CELL_REAR = CELL_REAR;
    CELL_output.TYPE=TYPE;
    CELL_output.SUBMOD_IND = Submod_ind;
    %---
    cd(current_path);
    disp('Wavelength and angle dependent absoptance calculated.')
end
end