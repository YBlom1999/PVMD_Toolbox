function [HeatGen,ind_Diode_new] = CalculateGeneratedHeat(TOOLBOX_input,I,Impp,IVParameters,ind_Diode,Ncells,days,T)
%CalculateGeneratedHeat calculates the generated heat due to reverse bias
%operating point
%
% This function calculates how much heat should be added to the thermal
% calculation
% 
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% I : double
%   The different current for the IV curves
% Impp : double
%   The maximum powerpoint current througout the year
% Parameters : double
%   Parameters of the equivalent circuit
% ind_Diode : logical
%   The index of which bypass diodes are active
% Ncells : double
%   The number of cells
% days : double
%   Index on which days irradiance is received
% T : double
%   The temperature of all cells
%
% Returns
% -------
% HeatGen : double
%   The additional heat that is generated
% ind_Diode_new: logical
%   An updated index of which bypass diodes are active
%
% Developed by by Youri Blom

HeatGen = zeros(length(Impp),Ncells);
Nby = size(ind_Diode,1);
K=1.3806e-23;
q=1.6022e-19;

ind_Diode_new = zeros(Nby,length(Impp));

for t= 1:length(days)
    Impp_t = Impp(days(t));
    if Impp_t==0; continue; end
    ind_I = I == Impp_t;
    Nby_act = ind_Diode(:,days(t),ind_I);
    ind_Diode_new(:,days(t)) = Nby_act;
    if strcmp(TOOLBOX_input.electric.ModuleType,'Butterfly')
        param_ind = (1:Ncells)+(t-1)*Ncells;
        Impp_t_full = CalculateCurrentButterfly(TOOLBOX_input,I,Impp_t,IVParameters(param_ind,:),T(days(t),:));
    end
    for cell_i = 1:Ncells
        if Nby_act(ceil(Nby*cell_i/72)); continue; end
        param_ind = (t-1)*Ncells+cell_i;
        
        if sum(IVParameters(param_ind,1:5)) == 0; continue; end

        if strcmp(TOOLBOX_input.electric.ModuleType,'Butterfly')
            Impp_t = Impp_t_full(cell_i);
        end
        

        V = [-50:-1 0:5e-3:1.5];
        I_ = CalculateIVcurve(TOOLBOX_input,IVParameters(param_ind,:),T(days(t),cell_i),V);
        V_cell = interp1(I_,V,Impp_t,'linear','extrap');
        HeatGen(days(t),cell_i) = -V_cell*Impp_t;

        

    end
end

end

function Impp_full = CalculateCurrentButterfly(TOOLBOX_input,I,Impp,IVParameters,T)
%CalculateCurrentButterfly calculates the current in each cell at MPP for
%the butterfly topology
%
% This function calculates how much current flows through each cell, which
% is needed for the calculation of the heat generation
% 
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% I : double
%   The current
% Impp : double
%   The maximum powerpoint current
% IVParameters : double
%   Parameters of the equivalent circuit
% T : double
%   The temperature of all cells
%
% Returns
% -------
% Impp_full : double
%   The current for each cell at MPP
%
% Developed by by Youri Blom

Impp_full = zeros(length(T),1);
N_rows = TOOLBOX_input.Scene.module_mounting.CellRows;
V = [-50:-1 0:5e-3:1.5];
Nby = TOOLBOX_input.electric.numBypassDiodes;
for i = 1:Nby
    V_string1 = zeros(length(I),length(IVParameters)/2/Nby);
    V_string2 = zeros(length(I),length(IVParameters)/2/Nby);
    Start_ind = 2*(i-1)*N_rows; 
    ind1 = Start_ind + [1:0.5*N_rows, N_rows+1:1.5*N_rows];
    ind2 = Start_ind + [0.5*N_rows+1:N_rows, 1.5*N_rows+1:2*N_rows];

    for cell_i = 1:length(ind1)
        if sum(IVParameters(ind1(cell_i),1:5)) == 0
            IVParameters(ind1(cell_i),1) = 0;
            IVParameters(ind1(cell_i),2) = 1e-8;
            IVParameters(ind1(cell_i),3) = 1e8;
            IVParameters(ind1(cell_i),4) = 2;
            IVParameters(ind1(cell_i),5) = 1e-8;
        end
        I_ = CalculateIVcurve(TOOLBOX_input,IVParameters(ind1(cell_i),:),T(ind1(cell_i)),V);
        intp_ind = ~isinf(I_);
        V_string1(:,cell_i) = interp1(I_(intp_ind),V(intp_ind),I,'linear','extrap');

        if sum(IVParameters(ind2(cell_i),1:5)) == 0
            IVParameters(ind2(cell_i),1) = 0;
            IVParameters(ind2(cell_i),2) = 1e-8;
            IVParameters(ind2(cell_i),3) = 1e8;
            IVParameters(ind2(cell_i),4) = 2;
            IVParameters(ind2(cell_i),5) = 1e-8;
        end

        I_ = CalculateIVcurve(TOOLBOX_input,IVParameters(ind2(cell_i),:),T(ind2(cell_i)),V);
        intp_ind = ~isinf(I_);
        V_string2(:,cell_i) = interp1(I_(intp_ind),V(intp_ind),I,'linear','extrap');
        

    end
    V_string1 = sum(V_string1,2);
    V_string2 = sum(V_string2,2);
    
    V_steps = linspace(min([V_string1;V_string2]),max([V_string1;V_string2]),1000);
    I_string1 = interp1(V_string1,I,V_steps,"linear","extrap");
    I_string2 = interp1(V_string2,I,V_steps,"linear","extrap");

    I_total = I_string1 + I_string2;
    V_string_needed =  interp1(I_total,V_steps,Impp);

    if isfield(TOOLBOX_input.electric,'TaylorParam')
        [~,indmpp_string] = max(I_total.*V_steps);
        Vmpp_string = V_steps(indmpp_string);
        V_string_needed = max(V_string_needed,Vmpp_string);
    end
    Impp_full(ind1) = interp1(V_string1,I,V_string_needed,"linear","extrap");
    Impp_full(ind2) = interp1(V_string2,I,V_string_needed,"linear","extrap");
end


end


function I = CalculateIVcurve(TOOLBOX_input,IVParameters,T,V)
%CalculateIVcurve calculates the IV curve of a certain cell
%
% This function calculates the IV curve for a cell based on its circuit
% parameters
% 
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% IVparameters : double
%   The parameters of equivalent circuit
% T : double
%   The temperature
% V : double
%   The voltage axis
%
% Returns
% -------
% I : double
%   The current of the cell
%
% Developed by by Youri Blom
q = TOOLBOX_input.constants.q;
k = TOOLBOX_input.constants.k;
Vth=k*T/q;

if size(IVParameters,2) == 5
    reverseModel = 0;
elseif size(IVParameters,2) == 9
    reverseModel = 1;
elseif size(IVParameters,2) == 8
    reverseModel = 2;
end

Iph = IVParameters(1);
Rs = IVParameters(2);
Rsh = IVParameters(3);
n = IVParameters(4);
I0 = IVParameters(5);

if reverseModel == 1
    Be = IVParameters(6);
    phi_t = IVParameters(7);
    V_b = IVParameters(8);
    c = IVParameters(9);
elseif reverseModel == 2
    I0_rev = IVParameters(6);
    Vs_rev = IVParameters(7);
    Rsh_rev = IVParameters(8);
end

%Current at forward bias
z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+V)./(n*Vth*(1+Rs/Rsh)));
I=(Iph+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;

%Current at reverse bias (when a different model is considered)
ind_neg = V < 0;
if reverseModel == 1
    ind_inf = V < V_b;
    I_n = Iph-V_test/Rsh+c*V_test.^2;
    K_e = exp(Be*(1-sqrt((phi_t-V_b)./(phi_t-V_test))));
    I(ind_neg) = I_n./(1-K_e);
    I(ind_inf) = -1e5*V(ind_inf);
elseif reverseModel == 2
    I(ind_neg) = Iph+I0_rev*(exp((-V_test)/(Vs_rev))-1)-V_test/Rsh_rev;
end
end