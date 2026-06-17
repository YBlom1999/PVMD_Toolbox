function [TOOLBOX_input] = get_LCOE_residential_input(TOOLBOX_input)
% get_LCOE_residential_input obtains the input needed for the LCOE of the
% residential system calculations
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
%
% Returns
% -------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
%
% Author: Y. Blom


%Ask for heat pump and battery size
Q = {'Heat pump [kW]','Battery size [Wh]'};
default = {'0','0'};
Answers = str2double(inputdlg(Q,'Size system',1,default));
TOOLBOX_input.FinancialPart.Size_HP = Answers(1);
TOOLBOX_input.FinancialPart.Size_Bat = Answers(2);

current_path = pwd;
[~,~,data_folder] = get_folder_structure;
cd(fullfile(data_folder, 'Financial','NproData'));    %go to 'Types' folder (where cell types are stored)
list = dir('*.xlsx');                              %list the file names there
Q ='Choose Npro data';                    %ask user
cd(current_path);
A= listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name}); %user choice
TOOLBOX_input.FinancialPart.NproFile = list(A).name;

end