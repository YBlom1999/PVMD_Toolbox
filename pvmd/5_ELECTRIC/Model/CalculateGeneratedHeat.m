function [HeatGen,ind_Diode_new] = CalculateGeneratedHeat(I,Impp,Parameters,ind_Diode,Ncells,days,T,TOOLBOX_input)
%CalculateGeneratedHeat calculates the generated heat due to reverse bias
%operating point
%
% This function calculates how much heat should be added to the thermal
% calculation
% 
%
% Parameters
% ----------
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
% reverseParam : double
%   The parameters of the reverse breakdown model
% TOOLBOX_input : struct
%   Simulation parameters
%
% Returns
% -------
% HeatGen : double
%   The additional heat that is generated
%
% Developed by by Youri Blom

HeatGen = zeros(length(Impp),Ncells);
Nby = size(ind_Diode,1);

ind_Diode_new = zeros(Nby,length(Impp));

for t= 1:length(days)
    Impp_t = Impp(days(t));
    if Impp_t==0; continue; end
    [~,ind_I] = min(abs(I -Impp_t));
    Nby_act = ind_Diode(:,days(t),ind_I);
    ind_Diode_new(:,days(t)) = Nby_act;
    if strcmp(TOOLBOX_input.electric.TYPE,'Butterfly')
        param_ind = (1:Ncells)+(t-1)*Ncells;
        Impp_t_full = CalculateCurrentButterfly(I,Impp_t,Parameters(param_ind,:),T(days(t),:),TOOLBOX_input);
    end
    for cell_i = 1:Ncells
        if Nby_act(ceil(Nby*cell_i/Ncells)); continue; end
        param_ind = (t-1)*Ncells+cell_i;
        Iph = Parameters(param_ind,1);
        Rs = Parameters(param_ind,2);
        Rsh = Parameters(param_ind,3);
        n = Parameters(param_ind,4);
        I0 = Parameters(param_ind,5);

        if Iph+Rs+Rsh+n+I0 == 0; continue; end

        if strcmp(TOOLBOX_input.electric.TYPE,'Butterfly')
            Impp_t = Impp_t_full(cell_i);
        end

        

        V = [-50:-1 0:5e-3:1.5];
        I_ = CalculateIVcurve(Iph,Rs,Rsh,n,I0,T(days(t),cell_i),TOOLBOX_input.electric.reverseParam,V);
        V_cell = interp1(I_,V,Impp_t,'linear','extrap');
        HeatGen(days(t),cell_i) = -V_cell*Impp_t;
        

    end
end

end


function Impp_full = CalculateCurrentButterfly(I,Impp,Parameters,T,TOOLBOX_input)
%CalculateCurrentButterfly calculates the current in each cell at MPP for
%the butterfly topology
%
% This function calculates how much current flows through each cell, which
% is needed for the calculation of the heat generation
% 
%
% Parameters
% ----------
% I : double
%   The current
% Impp : double
%   The maximum powerpoint current
% Parameters : double
%   Parameters of the equivalent circuit
% T : double
%   The temperature of all cells
% TOOLBOX_input : struct
%   Simulation parameters
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
reverseParam = TOOLBOX_input.electric.reverseParam;
for i = 1:Nby
    V_string1 = zeros(length(I),length(Parameters)/2/Nby);
    V_string2 = zeros(length(I),length(Parameters)/2/Nby);
    Start_ind = 2*(i-1)*N_rows; 
    ind1 = Start_ind + [1:0.5*N_rows, N_rows+1:1.5*N_rows];
    ind2 = Start_ind + [0.5*N_rows+1:N_rows, 1.5*N_rows+1:2*N_rows];

    for cell_i = 1:length(ind1)
        Iph = Parameters(ind1(cell_i),1);
        Rs = Parameters(ind1(cell_i),2);
        Rsh = Parameters(ind1(cell_i),3);
        n = Parameters(ind1(cell_i),4);
        I0 = Parameters(ind1(cell_i),5);
        if Iph+Rs+Rsh+n+I0 == 0
            Iph = 0;
            Rs = 1e-8;
            Rsh = 1e8;
            n = 2;
            I0 = 1e-8;

        end
        I_ = CalculateIVcurve(Iph,Rs,Rsh,n,I0,T(ind1(cell_i)),reverseParam,V);
        intp_ind = ~isinf(I_);
        V_string1(:,cell_i) = interp1(I_(intp_ind),V(intp_ind),I,'linear','extrap');

        Iph = Parameters(ind2(cell_i),1);
        Rs = Parameters(ind2(cell_i),2);
        Rsh = Parameters(ind2(cell_i),3);
        n = Parameters(ind2(cell_i),4);
        I0 = Parameters(ind2(cell_i),5);
        if Iph+Rs+Rsh+n+I0 == 0
            Iph = 0;
            Rs = 1e-8;
            Rsh = 1e8;
            n = 2;
            I0 = 1e-8;

        end
        I_ = CalculateIVcurve(Iph,Rs,Rsh,n,I0,T(ind2(cell_i)),reverseParam,V);
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


function I = CalculateIVcurve(Iph,Rs,Rsh,n,I0,T,reverseParam,V)
%CalculateIVcurve calculates the IV curve of a certain cell
%
% This function calculates the IV curve for a cell based on its circuit
% parameters
% 
%
% Parameters
% ----------
% Iph : double
%   The photo generated current
% Rs : double
%   The series resistance
% Rsh : double
%   The shunt resistance current
% n : double
%   The ideality factor
% I0 : double
%   The saturation current
% T : double
%   The temperature
% reverseParam : double
%   The parameters of the reverse breakdown model
% V : double
%   The voltage axis
%
% Returns
% -------
% Impp_full : double
%   The current for each cell at MPP
%
% Developed by by Youri Blom
K=1.3806e-23;
q=1.6022e-19;
Vth=K*T/q;

z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+V)./(n*Vth*(1+Rs/Rsh)));
I=(Iph+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;

% Apply different model for negative currents
if ~isnan(reverseParam)
    Be = reverseParam(1);
    phi_t = reverseParam(2);
    V_b = reverseParam(3);
    c = reverseParam(4);
    ind_neg = V < 0;
    ind_inf = V < V_b;
    I_n = Iph-V(ind_neg)/Rsh+c*V(ind_neg).^2;
    K_e = exp(Be*(1-sqrt((phi_t-V_b)./(phi_t-V(ind_neg)))));
    I(ind_neg) = I_n./(1-K_e);
    I(ind_inf) = -1e5*V(ind_inf);
end

end