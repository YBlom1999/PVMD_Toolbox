function [V_,day,night,param] = Simulated_IV(Acell,J_,T_,Op,numCells,I,IVtype,InputElectric,degType,Losses) %changed by youri
%%Generated current to over cell IV curves
%Author: Abdallah Nour El Din, Malte Vogt(Reworked)
%Needs to be commented
%Input
%Acell cell area in [m2]
%J_-Implied photocurrent in every every absorber layer in number of electrons per m2 format time in rows
%and cells in collumns
%T_-Temperature for each cell time in rows and cells in collumns
%Op - Optical shading loss in
%numCells-number of cells
%IVtype - the file name of a file from the data-folder, the file needs to
%contain the Temperature and Irradiance dependence of the IV curve parameters

%% Running current mapping
tic
%Reduce current based on shading
J_=J_*(1-Op*0.01);

%Sum over all cells
j=sum(J_,2);
%timesteps with zero current (night)
night=find(j==0);
%timesteps with current greater than zero (day)
day=find(j>0);
%delete timesteps (rows) with zero current
J_(night,:)=[];
T_(night,:)=[];
K=1.3806e-23;
q=1.6022e-19;

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
[src_folder,~,~] = get_folder_structure;
parameter_file = fullfile(src_folder,'5_ELECTRIC','Model','Data',IVtype);
load(parameter_file, 'J_sc_J','J_sc_T','J0_J','J0_T','R_s_J','R_s_T','R_sh_J','R_sh_T','n_J','n_T');


%IV curves for each condition with resolution of I
sim=zeros(length(cond(:,1)),length(I));
%Voltage steps
V=[-50:-1 0:5e-3:1.5];

%exporting diode parameters for testing purposes
export_param=zeros(length(cond(:,1)),5);

%Calculate IV-curves for all bins
for i=1:length(cond(:,1))
    J=cond(i,1)*J_bin-0.5*J_bin;
    T=cond(i,2)*T_bin-0.5*T_bin+Tmin;

    %prepareing parameters for one diode model solved by lambert-wfunction
    Iph=Acell*polyval(J_sc_T,T).*polyval(J_sc_J,J)./polyval(J_sc_T,298.15);
    Vth=K*T/q;
    n=polyval(n_T,T).*polyval(n_J,J)./polyval(n_T,298.15);
    %Correction for cell area
    %due to fitting JV (instead of IV) from ASA
    %Rs & Rsh would be m2*Ohm and I0 would be in A/m2 without corrrection
    Rsh=(1/Acell)*exp(polyval(R_sh_T,T)).*exp(polyval(R_sh_J,J))./exp(polyval(R_sh_T,298.15));
    Rs=(1/Acell)*exp(polyval(R_s_T,T)).*exp(polyval(R_s_J,J))./exp(polyval(R_s_T,298.15));
    I0=Acell*exp(polyval(J0_T,T)).*exp(polyval(J0_J,J))./exp(polyval(J0_T,298.15));
    
    %adjust parameters with factors
    if degType == 1
        factors = InputElectric.factors_bot;
    elseif degType == 2
        Parameters = [Iph,Rs,Rsh,n,I0];
        [factors] = FactorsPerovskite(Parameters,InputElectric.k_needed_top,InputElectric.Scenario,T);
    end
    Iph = Iph*factors(1);
    Rs = Rs*factors(2);
    Rsh = Rsh*factors(3);
    n = n*factors(4);
    I0 = 10^(log10(I0)*factors(5));

    %one diode model solved by lambert-wfunction
    z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+V)./(n*Vth*(1+Rs/Rsh)));
    I_=(Iph+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;

    % Apply different model for negative currents
    if ~isnan(InputElectric.reverseParam)
        Be = InputElectric.reverseParam(1);
        phi_t = InputElectric.reverseParam(2);
        V_b = InputElectric.reverseParam(3);
        c = InputElectric.reverseParam(4);
        ind_neg = V < 0;
        ind_inf = V < V_b;
        I_n = Iph-V(ind_neg)/Rsh+c*V(ind_neg).^2;
        K_e = exp(Be*(1-sqrt((phi_t-V_b)./(phi_t-V(ind_neg)))));
        I_(ind_neg) = I_n./(1-K_e);
        I_(ind_inf) = -1e5*V(ind_inf);
    end

    %Interpolation to comon currents
    ind = isinf(I_).*I_<0;
    I_(ind) = -1e5*V(ind);
    sim(i,:)=interp1(I_,V,I,'linear','extrap');

    %exporting diode parameters for testing purposes
    export_param(i,:)=[Iph,Rs,Rsh,n,I0];
end
%exporting diode parameters for testing purposes
% save('diodeParameter.mat','export_param');

%% Remap the bins in the time slots
V_=zeros(length(day)*numCells,length(I));
param = zeros(length(day)*numCells,5); %added by Youri
for i=1:length(cond(:,1))
    index=find(pointer(:,1)==cond(i,1) & pointer(:,2)==cond(i,2));
    V_(index,:)=repmat(sim(i,:),length(index),1);
    param(index,:) = repmat(export_param(i,:),length(index),1); %added by Youri
end

%% Update Performance after loss
if isstruct(Losses)
    [param,V_] = AccountLosses(Losses,V_,I,param,T_,V,InputElectric);

end

end

function [param_new,V_new] = AccountLosses(Losses,V,I,param,T,V_axis,InputElectric)
% AccountLosses accounts for the losses that the cell experience due to
% degradation
%
% It uses the calculated loss values and adjust the IV curves accordingly
%
% Parameters
% ----------
% Losses: struct
%   The calculated values for the different losses
% V: double
%   The voltage of all IV curves
% I: double
%   The current of all IV curves
% param: double
%   The equivalent circuit parameters of all cells
% T: double
%   The temperature of all cells
%
% Returns
% -------
% param_new: double
%   The new values for the equivalent circuit parameters
% V_new: double
%   The new values for the voltage of the IV curves

cell_index = repmat(1:144,1,size(param,1)/144);
param_new = zeros(size(param));
V_new = zeros(size(V));

K=1.3806e-23;
q=1.6022e-19;
Vth=K*T/q;

for i = 1:length(cell_index)
    Iph = param(i,1);
    Rs = param(i,2);
    Rsh = param(i,3);
    n = param(i,4);
    I0 = param(i,5);

    if Iph == 0; continue; end

    k_Isc = (1-Losses.k_dis(cell_index(i)));
    k_Voc = (1-Losses.k_mois(cell_index(i)));
    k_Rs = Losses.k_TC(cell_index(i));

    Iph_new = Iph*k_Isc;
    I0_new = (I0^k_Voc)*Iph^(1-k_Voc);

    [Pmpp,mpp_ind] = max(V(i,:).*I);
    Impp = I(mpp_ind);
    
    if Impp> 0; Rs_new = Rs+Pmpp*k_Rs/Impp^2; else Rs_new =  Rs; end

    %one diode model solved by lambert-wfunction
    z=(Rs_new*I0_new/(n*Vth(i)*(1+Rs_new/Rsh)))*exp((Rs_new*(Iph_new+I0_new)+V_axis)./(n*Vth(i)*(1+Rs_new/Rsh)));
    I_=(Iph_new+I0_new-V_axis/(Rsh))/(1+Rs_new/Rsh)-lambertw(z).*(n*Vth(i))/Rs_new;

    if ~isnan(InputElectric.reverseParam)
        Be = InputElectric.reverseParam(1);
        phi_t = InputElectric.reverseParam(2);
        V_b = InputElectric.reverseParam(3);
        c = InputElectric.reverseParam(4);
        ind_neg = V_axis < 0;
        ind_inf = V_axis < V_b;
        I_n = Iph-V_axis(ind_neg)/Rsh+c*V_axis(ind_neg).^2;
        K_e = exp(Be*(1-sqrt((phi_t-V_b)./(phi_t-V_axis(ind_neg)))));
        I_(ind_neg) = I_n./(1-K_e);
        I_(ind_inf) = -1e5*V_axis(ind_inf);
    end
    ind_inf = ~isinf(I_);
    V_new(i,:) = interp1(I_(ind_inf),V_axis(ind_inf),I,'linear','extrap');
    param_new(i,:) = [Iph_new,Rs_new,Rsh,n,I0_new];
end

end