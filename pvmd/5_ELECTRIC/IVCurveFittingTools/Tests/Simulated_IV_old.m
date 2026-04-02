function [V_,day,night] = Simulated_IV_old(Acell,J_,T_,Op,numCells,I,IVtype)
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
J_bin=0.4;
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
IVfilePath=join(['Data\', IVtype]);
load(IVfilePath);

%IV curves for each condition with resolution of I
sim=zeros(length(cond(:,1)),length(I));
%Voltage steps
V=[-15:-1 0:5e-3:1.5];

%exporting diode parameters for testing purposes
%export_param=zeros(length(cond(:,1)),5);

%Calculate IV-curves for all bins
for i=1:length(cond(:,1))
        J=cond(i,1)*J_bin-0.5*J_bin; 
        T=cond(i,2)*T_bin-0.5*T_bin+Tmin; 
        %Temperature below -30deg (243K) can cause errors thus manual limt (possible cause might be that ASA data
        %used in fit does not reach that far)
        %if T<243
         %   disp('Warning! Temperature below -30 degC (243K) this can cause errrors.');
         %   disp('For purposes of electrical simulation temperature is thus correct to -30degC (243K) for those intances.')
         %   T=243;
        %end  
        %prepareing parameters for one diode model solved by lambert-wfunction
        Iph=Acell*J;%.*(1+3e-4.*(T-298.15)); %Increased absorption with temperature only for Si not perovskite
        Vth=K*T/q;
        n=polyval(n_T,T).*polyval(n_J,J)./polyval(n_T,298.15);
        %Correction for cell area
        %due to fitting JV (instead of IV) from ASA
        %Rs & Rsh would be m2*Ohm and I0 would be in A/m2 without corrrection  
        Rsh=(1/Acell)*exp(polyval(R_sh_T,T)).*exp(polyval(R_sh_J,J))./exp(polyval(R_sh_T,298.15));
        Rs=(1/Acell)*exp(polyval(R_s_T,T)).*exp(polyval(R_s_J,J))./exp(polyval(R_s_T,298.15));
        I0=Acell*exp(polyval(J0_T,T)).*exp(polyval(J0_J,J))./exp(polyval(J0_T,298.15)) ;
        %one diode model solved by lambert-wfunction
        z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+V)./(n*Vth*(1+Rs/Rsh)));
        I_=(Iph+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;
        %Interpolation to comon currents
        sim(i,:)=interp1(I_,V,I);
        
        %exporting diode parameters for testing purposes
        %export_param(i,:)=[Iph,Rs,Rsh,n,I0];
end
%exporting diode parameters for testing purposes
%save diodeParameter.mat export_param;

%% Remap the bins in the time slots
V_=zeros(length(day)*numCells,length(I));
for i=1:length(cond(:,1))
    index=find(pointer(:,1)==cond(i,1) & pointer(:,2)==cond(i,2));
    V_(index,:)=repmat(sim(i,:),length(index),1);
end

end

