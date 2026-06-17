function [k_pred,CELL_output,MODULE_output,WEATHER_output,THERMAL_output] = DegDiscoloration(TOOLBOX_input,T,UV,A,Ea)
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
%
% Returns
% -------
% k_pred: double
%   The predicted degradation rate due to discoloration
% CELL_output : struct
%   Simulation results of the CELL module
% MODULE_output : struct
%   Simulation results of the MODULE module
% WEATHER_output : struct
%   Simulation results of the WEATHER module
% THERMAL_output : struct
%   Simulation results of the THERMAL module
%
% Developed by by Youri Blom
k = TOOLBOX_input.constants.k;
q = TOOLBOX_input.constants.q;

%load degradation file
DiscolorationFilePath = [TOOLBOX_input.Degradation.Discoloration_file,'.xlsx'];
DiscolorationData = readcell(DiscolorationFilePath,"Sheet","OpticalLoss","Range","A2:B8");
factor = cell2mat(DiscolorationData(:,1));
Jsc_loss = cell2mat(DiscolorationData(:,2));

k_pred = A*UV.*exp(-Ea*q/k./T);

factor_needed = interp1(Jsc_loss,factor,sum(k_pred),"linear","extrap");
[CELL_output,TOOLBOX_input] = update_CELL_sim(TOOLBOX_input,factor_needed);
[MODULE_output] = MODULE_main(TOOLBOX_input, CELL_output);
[WEATHER_output] = WEATHER_main(TOOLBOX_input, MODULE_output);
[THERMAL_output] = THERMAL_main(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output);


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
DiscolorationFilePath = [TOOLBOX_input.Degradation.Discoloration_file,'.xlsx'];
DiscolorationData = readcell(DiscolorationFilePath,"Sheet","Settings","Range","B2:B3");
nk_name = 'Polyolefin-UVT';%DiscolorationData{1,1};
loc_nkfile = which([nk_name,'.xlsx']);
nk_name_new = strrep(loc_nkfile,nk_name,'DiscEncapsulation');
d = DiscolorationData{2,1};

DiscolorationData = readcell(DiscolorationFilePath,"Sheet","ChangeEQE","Range","A2:D182");
wav = cell2mat(DiscolorationData(:,1));
ParAbs_init = cell2mat(DiscolorationData(:,2));
ParAbs_final = cell2mat(DiscolorationData(:,3));
R_final = cell2mat(DiscolorationData(:,4));

Ni=read_nk(nk_name,wav*1e3);
n_orig = real(Ni);
k_orig = imag(Ni);

Delta_ParAbs = factor_needed*(ParAbs_final - ParAbs_init);

Delta_alpha = -log((1-Delta_ParAbs-R_final)./(1-R_final))/d;
Delta_k_EVA = max(Delta_alpha.*wav*1e-9/4/pi,0);

nm = wav;
n = n_orig'+factor_needed;
k = k_orig'+Delta_k_EVA;

text = table(nm,n,k);
writetable(text,nk_name_new,"Sheet","nk","Range","A1:C182")
addpath(erase(nk_name_new,'DiscEncapsulation.xlsx'))

%Update input file
GP_input_file = which([TOOLBOX_input.deviceOptic.GenProFile,'.xlsx']);
new_GP_input_file = 'New_GP_input.xlsx';

copyfile(GP_input_file,new_GP_input_file);

cellStructure = readtable(new_GP_input_file,"Sheet","Layers");

line_i = 1;
while true
    layer_i = cellStructure{line_i,1};
    if strcmp(layer_i,nk_name)
        cellStructure{line_i,1} = {'DiscEncapsulation'};
        writetable(cellStructure,new_GP_input_file,"Sheet","Layers");
        break
    end
    if line_i > 100
        warning('Layer not found')
        break
    end
    line_i = line_i+1;
end


% Redo Optical simulation
TOOLBOX_input.deviceOptic.GenProFile = erase(new_GP_input_file,'.xlsx');
TOOLBOX_input.script = 1;
TOOLBOX_input.electric.electricplot = 0;
CELL_output = CELL_main(TOOLBOX_input);

% Delete files
delete(new_GP_input_file);
% delete(nk_name_new);
end