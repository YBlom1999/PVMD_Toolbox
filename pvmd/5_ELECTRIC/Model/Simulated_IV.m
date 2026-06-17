function [Vcell,day,night,param] = Simulated_IV(TOOLBOX_input,Acell,numCells,Shading,IVtype,J_full,T_full,Irr,I,cell_ind,degType)
% Simulated_IV calculates the cell IV curves with the calibrated lumped
% element method for both forward and reverse bias
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% Acell: double
%   The area of the cell
% numCells : double
%   The number of cells
% Shading: double
%   The optical shading due to the contacts
% IVtype: string
%   The IV characteristics of the cell
% J_full: double
%   The photoreceived of all cells
% T_full: double
%   The temperature of all cells
% Irr : double
%   The received irradiance of all cells
% cell_ind : double
%   The index of which cell is simulated
% degType: double
%   An indicator of how the degradation should be simulated
%
% Returns
% -------
% Vcell : double
%   The voltage of the cell IV curves
% day: double
%   The moments at which there is daylight
% night: double
%   The moments at which it is night
% param : double
%   The parameters of the equivalent circuit
%
% Written by Abdallah Nour El Din
% Adjusted by M. R. Vogt & Y. Blom

% Read the settings
k_mois = TOOLBOX_input.electric.deg_sil.k_mois;
factor_LID = TOOLBOX_input.electric.deg_sil.factor_LID;
filename_deg = TOOLBOX_input.electric.deg_sil.filename;
k_needed = TOOLBOX_input.electric.k_needed;
degScenario = TOOLBOX_input.electric.degScenario;

%Reduce current based on shading
J_full=J_full*(1-Shading*0.01);

%Sum over all cells
j=sum(J_full,2);
%timesteps with zero current (night)
night=find(j==0);
%timesteps with current greater than zero (day)
day=find(j>0);
%delete timesteps (rows) with zero current
J_ = J_full; J_(night,:)=[];
T_ = T_full; T_(night,:)=[];

q = TOOLBOX_input.constants.q;
k = TOOLBOX_input.constants.k;

Tmin=min(min(T_));

%Turn matrix into a vector
J_=reshape(J_',numel(J_),1);
T_=reshape(T_',numel(T_),1);
%Mapp cells into 0.4 A/m2 and 0.3 degC bins
J_bin=1;%0.4;
T_bin=0.3;
pointer=zeros(length(J_),2);
pointer(:,1)=ceil((J_-J_bin)/J_bin);
pointer(:,2)=ceil((T_-Tmin)/T_bin);
%Identify unique conditions bins
cond=unique(pointer,'rows');
%Removes conditions with J_ < J_bin A/m2
index=find(cond(:,1)>1,1);
cond=cond(index:end,:);

%% Caculating IV curves

%Load IV curves from ASA
[paramForwardBias,paramReverseBias,paramLightSoaking] = readIVcurveFile(IVtype);
Jph_J = paramForwardBias.Jph_J;
Jph_T = paramForwardBias.Jph_T;
J0_J = paramForwardBias.J0_J;
J0_T = paramForwardBias.J0_T;
n_J = paramForwardBias.n_J;
n_T = paramForwardBias.n_T;
Rsh_J = paramForwardBias.Rsh_J;
Rsh_T = paramForwardBias.Rsh_T;
Rs_J = paramForwardBias.Rs_J;
Rs_T = paramForwardBias.Rs_T;

%IV curves for each condition with resolution of I
sim=zeros(length(cond(:,1)),length(I));
%Voltage steps
V=[-15:-1 0:5e-3:1.5];

%exporting diode parameters for testing purposes
if paramReverseBias.modelType ~= 2
    export_param=zeros(length(cond(:,1)),5);
else
    export_param=zeros(length(cond(:,1)),8);
end

%Calculate IV-curves for all bins
for i=1:length(cond(:,1))
    J=cond(i,1)*J_bin-0.5*J_bin;
    T=cond(i,2)*T_bin-0.5*T_bin+Tmin;
    Vth=k*T/q;

    %prepareing parameters for one diode model solved by lambert-wfunction
    Iph=Acell*polyval(Jph_T,T).*polyval(Jph_J,J)./polyval(Jph_T,298.15);
    n=polyval(n_T,T).*polyval(n_J,J)./polyval(n_T,298.15);
    Rsh=(1/Acell)*exp(polyval(Rsh_T,T)).*exp(polyval(Rsh_J,J))./exp(polyval(Rsh_T,298.15));
    Rs=(1/Acell)*exp(polyval(Rs_T,T)).*exp(polyval(Rs_J,J))./exp(polyval(Rs_T,298.15));
    I0=Acell*exp(polyval(J0_T,T)).*exp(polyval(J0_J,J))./exp(polyval(J0_T,298.15));
    
    %adjust parameters with factors
    if degType == 1
        Parameters = [Iph,Rs,Rsh,n,I0];
        if k_mois == 0
            factors = [1,1,1,1,1];
        else
            [factors] = FactorsSilicon(Parameters,k_mois,T,filename_deg);
        end
        Iph = Iph*factors(1);
        Rs = Rs*factors(2);
        Rsh = Rsh*factors(3);
        n = n*factors(4);
        I0 = 10^(log10(I0)*factors(5))/factor_LID;
    elseif degType == 2
        Parameters = [Iph,Rs,Rsh,n,I0];
        [factors] = FactorsPerovskite(Parameters,k_needed{cell_ind},degScenario{cell_ind},T);
        Iph = Iph*factors(1);
        Rs = Rs*factors(2);
        Rsh = Rsh*factors(3);
        n = n*factors(4);
        I0 = 10^(log10(I0)*factors(5));
    end

    %one diode model solved by lambert-wfunction
    z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+V)./(n*Vth*(1+Rs/Rsh)));
    I_=(Iph+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw_pvmd(z).*(n*Vth)/Rs;

    % Apply different model for negative currents
    if paramReverseBias.modelType == 1
        Be = paramReverseBias.Be;
        phi_t = paramReverseBias.phi_t;
        V_b = paramReverseBias.V_b;
        c = paramReverseBias.c;
        ind_neg = V < 0;
        I_n = Iph-V(ind_neg)/Rsh+c*V(ind_neg).^2;
        K_e = exp(Be*(1-sqrt((phi_t-V_b)./(phi_t-V(ind_neg)))));
        I_(ind_neg) = I_n./(1-K_e);
    elseif paramReverseBias.modelType == 2
        J0_rev_J = paramReverseBias.J0_rev_J;
        J0_rev_T = paramReverseBias.J0_rev_T;
        Vs_rev_J = paramReverseBias.Vs_rev_J;
        Vs_rev_T = paramReverseBias.Vs_rev_T;
        Rsh_rev_J = paramReverseBias.Rsh_rev_J;
        Rsh_rev_T = paramReverseBias.Rsh_rev_T;

        Vs_rev=polyval(Vs_rev_T,T).*polyval(Vs_rev_J,J)./polyval(Vs_rev_T,298.15);
        Rsh_rev=(1/Acell)*exp(polyval(Rsh_rev_T,T)).*exp(polyval(Rsh_rev_J,J))./exp(polyval(Rsh_rev_T,298.15));
        I0_rev=Acell*exp(polyval(J0_rev_T,T)).*exp(polyval(J0_rev_J,J))./exp(polyval(J0_rev_T,298.15));
        
        ind_neg = V < 0;
        I_(ind_neg) = Iph+I0_rev*(exp((-V(ind_neg))/(Vs_rev))-1)-V(ind_neg)/Rsh_rev;
    end

    %Interpolation to comon currents
    ind = isinf(I_).*I_<0;
    I_(ind) = -1e5*V(ind);
    sim(i,:)=interp1(I_,V,I,'linear','extrap');

    if paramReverseBias.modelType == 0
        export_param(i,:)=[Iph,Rs,Rsh,n,I0];
    elseif paramReverseBias.modelType == 1
        export_param(i,:)=[Iph,Rs,Rsh,n,I0,Be,phi_t,V_b,c];
    elseif paramReverseBias.modelType == 2
        export_param(i,:)=[Iph,Rs,Rsh,n,I0,I0_rev,Vs_rev,Rsh_rev];
    end
end


%% Remap the bins in the time slots
Vcell=zeros(length(day)*numCells,length(I));
param = zeros(length(day)*numCells,size(export_param,2));
for i=1:length(cond(:,1))
    index=find(pointer(:,1)==cond(i,1) & pointer(:,2)==cond(i,2));
    Vcell(index,:)=repmat(sim(i,:),length(index),1);
    param(index,:) = repmat(export_param(i,:),length(index),1);
end

%% Account for lightsoaking
if paramLightSoaking.A ~= 0 && numel(Irr) > numCells
    A = paramLightSoaking.A;
    B = paramLightSoaking.B;
    Ea = paramLightSoaking.Ea;

    Irr_mean = mean(Irr,2);
    T_mean = mean(T_full,2);
    if length(Irr_mean) > 24
        Irr_day = reshape(Irr_mean,24,length(Irr_mean)/24)';
        T_day = reshape(T_mean,24,length(T_mean)/24)';
    else
        Irr_day = Irr_mean;
        T_day = T_mean;
    end

    Irr_cum = cumsum(Irr_day,2);
    tau = 1./(B*exp(-q*Ea/k./T_day));
    delta_Vmpp_day = A*exp(-Irr_cum./tau);

    delta_Vmpp = reshape(delta_Vmpp_day',length(Irr_mean),1);
    delta_Vmpp = delta_Vmpp(day);
    delta_Vmpp = repelem(delta_Vmpp,numCells);

    if size(Vcell,1) == size(delta_Vmpp,1)
        Vcell = Vcell - delta_Vmpp;
    else
        Vcell = Vcell - delta_Vmpp';
    end
end

end

