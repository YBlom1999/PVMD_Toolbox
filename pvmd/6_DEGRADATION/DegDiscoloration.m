function [k_pred,CELL_output,MODULE_output,WEATHER_output,THERMAL_output] = DegDiscoloration(TOOLBOX_input,T,UV,A,Ea,CONSTANTS)
%DegDiscoloration Calculates the degradation due to discoloration
%
% This function calculates the degradation rate caused by discoloration of
% the encapsulation
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
% T : double
%   Temperature of the module
% UV : double
%   The UV light received by the module
% A : double
%   The pre-exponential constant of discoloration
% Ea : double
%   Activation energy of the discoloration process
% CONSTANTS : struct
%   Structure of physical constants
%
% Returns
% -------
% k_pred: double
%   The predicted degradation rate due to discoloration
%
% Developed by by Youri Blom
k_b = CONSTANTS.k_b;
q = CONSTANTS.q;

%load degradation file
deg_file = TOOLBOX_input.Degradation.Discoloration_file;
[~,~,data_folder] = get_folder_structure;
filename =fullfile(data_folder,'Degradation','Discoloration',deg_file);
load(filename,'Jsc_loss','factor');


k_pred = A*UV.*exp(-Ea*q/k_b./T);

factor_needed = interp1(Jsc_loss,factor,sum(k_pred),"linear","extrap");
[CELL_output,TOOLBOX_input] = update_CELL_sim(TOOLBOX_input,factor_needed);
[MODULE_output, TOOLBOX_input] = MODULE_main(TOOLBOX_input, CELL_output);
[WEATHER_output, TOOLBOX_input] = WEATHER_main(TOOLBOX_input, MODULE_output);
[THERMAL_output, ~] = THERMAL_main(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output);


end


function [CELL_output,TOOLBOX_input] = update_CELL_sim(TOOLBOX_input,factor_needed)
% update_nk_file Updates the nk data of the encapsulant
%
% It changes the nk data of the encapsulant and stores it in the GenPro
% folder
%
% Parameters
% ----------
%
% Returns
% -------
% CELL_output : struct
%   The Cell simulation with the updated EVA
% TOOLBOX_input : struct
%   User input parameters to the toolbox
%
% Author: Youri Blom

% Update NK
deg_file = TOOLBOX_input.Degradation.Discoloration_file;
[pvmd_folder,~,data_folder] = get_folder_structure;
filename =fullfile(data_folder,'Degradation','Discoloration',deg_file);
load(filename,'d','ParAbs_final','R_final','ParAbs_init','wav','nk_name');

nk_file = fullfile(pvmd_folder,'1_CELL','GenPro4','nk',nk_name);
new_nk_file = fullfile(pvmd_folder,'1_CELL','GenPro4','nk','NK_updated.nk');
data = importdata(nk_file);
wav_orig = data.data(:,1);
n_orig = data.data(:,2);
k_orig = data.data(:,3);

Delta_ParAbs = factor_needed*(ParAbs_final - ParAbs_init);

Delta_alpha = -log((1-Delta_ParAbs-R_final)./(1-R_final))/d;
Delta_k_EVA = max(Delta_alpha.*wav*1e-9/4/pi,0);
[~,end_ind] = min(abs(wav_orig-wav(end)));
nm = wav_orig;
n = n_orig+factor_needed;
k = k_orig+[interp1(wav,Delta_k_EVA,wav_orig(1:end_ind),'linear','extrap');zeros(length(wav_orig)-end_ind,1)];

text = table(nm,n,k);
writetable(text,new_nk_file,'FileType','text')

%Update input file
GP_input_file = TOOLBOX_input.deviceOptic.GenProFile;
new_GP_input_file = 'New_GP_input.m';

current_path = pwd;
cd(fullfile(pvmd_folder,'1_CELL','Types'))
copyfile(GP_input_file,new_GP_input_file);


fileID_orig = fopen(GP_input_file,'r');
fileID = fopen(new_GP_input_file,'w');
line_i = 1;
while true
    line = fgetl(fileID_orig);
    if strcmp(line,'end')
        break
    end
    if line_i == 1
        line = strrep(line,erase(GP_input_file,'.m'),erase(new_GP_input_file,'.m'));
    end
    if contains(line,erase(nk_name,'.nk'))
        line = strrep(line,erase(nk_name,'.nk'),'NK_updated');
        
    end
    fprintf(fileID,'%4s\n',line);
    line_i = line_i+1;
end
fclose(fileID); fclose(fileID_orig);
cd(current_path);

% Redo Optical simulation
TOOLBOX_input.deviceOptic.GenProFile = new_GP_input_file;
TOOLBOX_input.script = 1;
TOOLBOX_input.electric.electricplot = 0;
[CELL_output,~] = CELL_main(TOOLBOX_input);

% Delete files
delete(fullfile(pvmd_folder,'1_CELL','Types',new_GP_input_file));
delete(new_nk_file)
end