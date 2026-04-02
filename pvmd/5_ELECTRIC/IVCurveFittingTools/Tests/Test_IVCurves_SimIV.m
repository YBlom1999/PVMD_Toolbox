%To test the IV curves in PVMD Toolbox
%Author Malte
clear;

%% set parameters for the test

%If true the Datasheet model is used else simulated IV part is called
TOOLBOX_input.electric.runDatasheet=false;

%specify parameters
if TOOLBOX_input.electric.runDatasheet
      %Applies to  single junction - Type defined in Cell parts determines, which
      %is used
      %Structure {'Voc [V]','Isc [A]','Vmpp [V]','Impp [A]','Kp [% / ?C]','Kv [% / ?C]','Ki [% / ?C]'};
      TOOLBOX_input.electric.datasheetValues = [40,8,35,7.5,-0.35,-0.30,0.03];
else
    %Specify a file from the 5_ELECTRIC\Model\data-folder, the file needs to
    %contain the Temperature and Irradiance dependence of the IV curve parameters   
    
    %Applies to  single junction modules
    TOOLBOX_input.electric.IVtype = 'Silicon_HZB2915_v4.mat';
end

%specify Loss due to shading by metalization[%]
TOOLBOX_input.electric.shading=0;
 

%% Set conditions 
numCells=1;
Acell=0.0246; %m2
T_ = (238.15:10:368.15); %First dimension temnperatures, second numCells
J_ = ones(length(T_),1)*358;%First dimension current densities in A/m2, second numCells

%% Account for Isc given in data sheet model
if TOOLBOX_input.electric.runDatasheet
    JscGenPro=10*max(evalin('base','CELL_output.CELL_FRONT.Jph'));%mA/cm2-->A/m2
    Conv_Datasheet_GenPro=TOOLBOX_input.electric.datasheetValues(2)/(JscGenPro*Acell); 
else
    Conv_Datasheet_GenPro=1;
end    

%maximum simulation current
Imax=max(J_,[],'all')*Acell*Conv_Datasheet_GenPro*1.05; 
%Resolution of the IV curve in A
I=0:5e-3:Imax;

[V_,day,night] = Simulated_IV(Acell,J_,T_, TOOLBOX_input.electric.shading,...
        numCells,I, TOOLBOX_input.electric.IVtype);
    
%% Plot IV curves
plot_IVcurve(V_,I,string(T_),"IV_Temp"+TOOLBOX_input.electric.IVtype);
plot_YvsTemp(T_,V_(:,1),'V_{oc} [V]',"Voc_Temp"+TOOLBOX_input.electric.IVtype);

%find Isc
ind_Isc=zeros(length(T_),1);
for h=1:length(T_)
    for i=1:length(I)-1
        if V_(h,i)>0&&V_(h,i+1)<0
            ind_Isc(h)=i;
            break
        end
    end
    if ind_Isc(h)==0
        [~,ind_Isc]=min(V_(h,:));
    end    
end
plot_YvsTemp(T_,I(ind_Isc),'I_{sc} [A]',"Isc_Temp"+TOOLBOX_input.electric.IVtype);

%find MPP
P_mpp=zeros(length(T_),1);
for h=1:length(T_)
    j=min(ind_Isc(h));
    I_forMPP=I(1:j);
    V_forMPP=V_(h,1:j);
    P_mpp(h)=max(V_forMPP.*I_forMPP,[],2);
end
plot_YvsTemp(T_,P_mpp,'P_{mpp} [W]',"Pmpp_Temp"+TOOLBOX_input.electric.IVtype);

%% check with ASA IV curves

%Define temperatures of the IV curves
T=T_-273.15;
numOfVariations=length(T);

%Line for starting and ending the import
startLine=7;
endLine=97;
numOfLines=endLine-startLine+1;

Voltage=zeros(numOfLines,numOfVariations);
Current=zeros(numOfLines,numOfVariations);

%Set file names
for i=1:numOfVariations
	filename=join(['JV_SHJ_HZB_2T_Bot_v7_',num2str(T(i)),'deg_1000Wm2.dat']);
    
   
    [Voltage(:,i), Current(:,i)] = importIVCurve(filename, [startLine, endLine]);
end

%source file gives negative voltage values
Voltage=-1*Voltage;

RMSE_J=zeros(numOfVariations,1);
NRMSE_J=zeros(numOfVariations,1);
RMSE_V=zeros(numOfVariations,1);
NRMSE_V=zeros(numOfVariations,1);
for i=1:numOfVariations
     [RMSE_J(i), NRMSE_J(i), RMSE_V(i),NRMSE_V(i) ] = ...
         plot_IVcurveDifference(V_(i,:),I/Acell,Voltage(:,i),Current(:,i),...
         ["SimIV "+T_(i),"ASA "+T_(i)],T_(i)+"k_IVDelta_"+TOOLBOX_input.electric.IVtype);
end
plot_YvsTemp(T_,NRMSE_J,'NRMSE J',"NRMSE_J_T_"+TOOLBOX_input.electric.IVtype);
plot_YvsTemp(T_,NRMSE_V,'NRMSE V',"NRMSE_V_T_"+TOOLBOX_input.electric.IVtype);
%check with values fitted to ASA
%load('diodeParameter.mat');
%plot(T_,export_param(:,1)./Acell);
%plot(T_,export_param(:,2).*Acell);
%plot(T_,export_param(:,3).*Acell);
%plot(T_,export_param(:,4));
%plot(T_,export_param(:,5)./Acell);


