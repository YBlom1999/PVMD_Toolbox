function [P_diffA, I_abs] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs,T_cell,V_opt,IVcurve_name)

%Properties of the module
N_cells = MODULE_output.N;
A_cell = MODULE_output.A;
Sf = TOOLBOX_input.electric.shading/100;
A_eff = N_cells*A_cell*(1-Sf);
%If the Isc is larger than the absorbed current (possible for SHJ), this
%needs to be taken from the parasitic absorption.
[src_folder,~,~] = get_folder_structure;
load(fullfile(src_folder,'5_ELECTRIC','Model','Data',IVcurve_name),'J_sc_J','J_sc_T');
I_sc=(A_eff/N_cells)*polyval(J_sc_T,T_cell).*polyval(J_sc_J,I_abs/(A_eff/N_cells))./polyval(J_sc_T,298.15);
Change_ParAbs = I_abs < I_sc;
P_diffA(Change_ParAbs) = P_diffA(Change_ParAbs)-N_cells*V_opt(Change_ParAbs).*(I_sc(Change_ParAbs)-I_abs(Change_ParAbs));
I_abs(Change_ParAbs) = I_sc(Chang