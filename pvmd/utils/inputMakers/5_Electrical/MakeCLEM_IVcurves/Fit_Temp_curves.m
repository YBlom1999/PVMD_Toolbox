function [J_sc_T,J_0_T,R_s_T,R_sh_T,n_T] = Fit_Temp_curves(Voltage, Current,T_levels)
%Fit_Temp_curves Calculates the dependencies of the equivalent circuit
%parameters w.r.t. the temperature
%
% This function reads the IV curves, fits it with the one-diode equivalent
% circuit and makes a fit w.r.t. the temperature
%
% Parameters
% ----------
% Voltage: double
%   The voltages of all the IV curves
% Current: double
%   The currents of all the IV curves
% T_levels : double
%   The temperature range for which IV curves are measured/simulated.
%
% Returns
% -------
% J_sc_T : double
%   The values of the fit for the short circuit current
% J_0_T : double
%   The values of the fit for the saturation current
% R_s_T : double
%   The values of the fit for the series resistance
% R_sh_T : double
%   The values of the fit for the shunt resistance
% n_T : double
%   The values of the fit for the ideality factor
%
% Developed by Malte Vogt and Youri Blom

N_curves = size(Voltage,2);
T_levels=T_levels+273.15;
par_T=zeros(N_curves,9);

color_begin = [1,0,0]; %Red
color_end = [0,0,1]; %Blue
color = [linspace(color_begin(1),color_end(1),N_curves)',linspace(color_begin(2),color_end(2),N_curves)',linspace(color_begin(3),color_end(3),N_curves)'];
figure
hold on; box on;

for j=1:N_curves
    V=Voltage(:,j);
    J=Current(:,j);
    
    
    %find Jsc
    ind_Jsc=find(~V);
    Jsc=J(ind_Jsc);
    
    %find Voc
    for i=1:length(V)-1
        if J(i)>0 && J(i+1)<0
            Voc=V(i)+((0-J(i))/(J(i+1)-J(i)))*(V(i+1)-V(i));
            break
        end
    end
    
    
    
    %cut off negative currents
    for i=1:length(J)-1
        if J(i)>0 && J(i+1)<0
            J=J(1:i+1);
            V=V(1:i+1);
            break
        end
    end
    
    %determine MPP
    P=V.*J;
    [M,i]=max(P);
    Vmpp=V(i);
    Jmpp=J(i);
    FF=M/(Voc*Jsc);
    
    
    % shunt
    f_rp=fit(V(ind_Jsc:ind_Jsc+5),J(ind_Jsc:ind_Jsc+5),'poly1');
    coeff_f_rp=coeffvalues(f_rp);
    Rp=-1.0/coeff_f_rp(1);
    
    
    rp=Rp/(Voc/Jsc);
    Vth=(1.3806*10^-23)*T_levels(j)/(1.602*10^-19);
    % x(1)=voc - x(2)=rs
    FF0=@(x) ((x(1) - log(x(1)+0.72))/(x(1)+1));
    FFs=@(x) FF0(x)*(1-1.1*x(2))+(x(2)^2)/5.4;
    FFsh=@(x) FFs(x)*(1-((x(1)+0.7)/x(1))*(FFs(x)/rp));
    fill=@(x) abs(FF - FF0(x));
    
    
    if j==1
        X0=Jsc/Voc;
    else
        X0=voc;
    end
    [voc,~]=fminsearch(fill,X0);
    %
    fill=@(x) abs(FF - FFsh(x));
    X0=[voc 1e-4];
    [voc,~]=fminsearch(fill,X0);
    n=Voc/(voc(1)*Vth);
    Rs=voc(2)*(Voc/Jsc);
    
    % saturation current
    D=[V J]; % x=J0
    for k=1:length(D(:,2))-1
        if D(k,2)<-1
            D=D(1:k+1,:);
            break
        end
    end
    
    Jmd=@(x,D) Jsc - x.*(exp((D(:,1)+D(:,2)*Rs)./(Vth*n))-1)-(D(:,1)+Rs.*D(:,2))./Rp;
    Err=@(x) sum((D(:,2)-Jmd(x,D)).^2);
    X0=1e-9;
    options=optimset('MaxIter',15000,'TolX',1e-9,'TolFun',1e-4);
    [J0,~]=fminsearch(Err,X0,options);
    
    
    par_T(j,:)=[Jsc, Voc, Rp, Rs, J0, n, FF, Jmpp, Vmpp];
    
    %%PLot
    Jmd=Jsc - J0.*(exp((V+J*Rs)./(Vth*n))-1)-(V+Rs.*J)./Rp;
    plot(V,J,'color',color(j,:))
    plot(V,Jmd,'--','color',color(j,:))
    
end
xlabel('Voltage [V]')
ylabel('Current density [A/m2]')
xlim([0, max(max(Voltage))])
ylim([0, 1.05*max(Current(1,:))]);

%% Fit shunt resistance
shunt=par_T(:,3)';
R_sh_T=polyfit(T_levels,log(shunt),4);

%plot
figure;
T_1K=min(T_levels):1:max(T_levels);
f=polyval(R_sh_T,T_1K);
f=exp(f);
hold on
plot(T_1K,f/1e2,'LineWidth',2)
plot(T_levels,shunt/1e2,'o','MarkerSize',8)
% grid on
%title('Shunt Resistance vs. Irradiance')
ylabel('Shunt Resistance [MOhomcm^2]')
xlabel('Temperature [K]')
hold off
box on
ax = gca;
ax.FontSize = 13;
%ylim([0 max(shunt)])
xlim([min(T_levels) max(T_levels)])
% RMSE
f1=exp(polyval(R_sh_T,T_levels));
Rel_mean_dif_Rsh=mean(abs(f1-shunt)./shunt);
NRMSE_Rsh= sqrt(sum((f1-shunt).^2)/(N_curves*mean(shunt)));
figuretext=sprintf('Fit of Rsh\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Rsh,Rel_mean_dif_Rsh);
text(300, f(1)/1e2, figuretext)


%% Fit series resitance
series=par_T(:,4)';
R_s_T=polyfit(T_levels,log(series),4);

%plot
figure;
f=polyval(R_s_T,T_1K);
f=exp(f);
hold on
plot(T_1K,f*1e4,'LineWidth',2);
plot(T_levels,series*1e4,'o','MarkerSize',8);
ylabel('Series Resistance [Ohm.cm^2]');
xlabel('Temperature [K]');
ax = gca;
ax.FontSize = 13;
xlim([min(T_levels) max(T_levels)])
%ylim([0 max(series)])
box on
hold off

% RMSE
f1=exp(polyval(R_s_T,T_levels));
Rel_mean_dif_Rs=mean(abs(f1-series)./series);
NRMSE_Rs= sqrt(sum((f1-series).^2)/(N_curves*mean(series)));
figuretext=sprintf('Fit of Rs\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Rs,Rel_mean_dif_Rs);
text(300, 0.8*f(1)*1e4, figuretext)

%% Fit J0
sat=par_T(:,5)';
J_0_T=polyfit(T_levels,log(sat),4);

%Plot
figure;
f=exp(polyval(J_0_T,T_1K));
hold on
plot(T_1K,f*1e8,'LineWidth',2);
plot(T_levels,sat*1e8,'o','MarkerSize',8);
xlabel('Temperature [K]');
ylabel('Saturation Current Density [pA/cm^2]');
ax = gca;
ax.FontSize = 13;
ax.YScale='log';
xlim([min(T_levels) max(T_levels)])
box on
hold off

% RMSE
f1=exp(polyval(J_0_T,T_levels));
Rel_mean_dif_J0=mean(abs(f1-sat)./sat);
NRMSE_J0= sqrt(sum((f1-sat).^2)/(N_curves*mean(sat)));
figuretext=sprintf('Fit of J0\nNRMSE=%.2f\nRMD=%.4f',NRMSE_J0,Rel_mean_dif_J0);
text(300, 0.8*f(1)*1e8, figuretext)


%% Fit idiality n
n=par_T(:,6)';
n_T=polyfit(T_levels,n,4);

% Plot
figure;
f=polyval(n_T,T_1K);
hold on
plot(T_1K,f,'LineWidth',2);
plot(T_levels,n,'o','MarkerSize',8);
hold off
ax = gca;
ax.FontSize = 13;
%ylim([0.8 1.2]);
xlim([min(T_levels) max(T_levels)]);
box on
xlabel('Temperature [K]');
ylabel('Ideality Factor');

%RMSE
f1=polyval(n_T,T_levels);
Rel_mean_dif_n=mean(abs(f1-n)./n);
NRMSE_n= sqrt(sum((f1-n).^2)/(N_curves*mean(n)));
figuretext=sprintf('Fit of idialty factor\nNRMSE=%.2f\nRMD=%.4f',NRMSE_n,Rel_mean_dif_n);
text(300, f(1), figuretext)

%% Fit fill factor
Fill=par_T(:,7)';
po=polyfit(T_levels,Fill,3);

%plot
figure;
hold on
f=polyval(po,T_1K);
plot(T_1K,f,T_levels,Fill,'o');
hold off

box on
xlabel('Temperature (K)');
ylabel('FF');

%RMSE
f1=polyval(po,T_levels);
Rel_mean_dif_Fill=mean(abs(f1-Fill)./Fill);
NRMSE_Fill= sqrt(sum((f1-Fill).^2)/(N_curves*mean(Fill)));
figuretext=sprintf('Fit of fill factor\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Fill,Rel_mean_dif_Fill);
text(300, 1.01*f(1), figuretext);

%% Fit Jsc
J_sc=par_T(:,1)';
J_sc_T=polyfit(T_levels,J_sc,4);

%plot
figure;
f=polyval(J_sc_T,T_1K);
hold on
plot(T_1K,f/10,'LineWidth',2);
plot(T_levels,J_sc/10,'o','MarkerSize',8);
ax = gca;
ax.FontSize = 13;
xlim([min(T_levels) max(T_levels)]);
hold off
box on
xlabel('Temperature (K)');
ylabel('Photo-generated Current Density [mA/cm^2]');

%Determine goodnes of fit
f1=polyval(J_sc_T,T_levels);
Rel_mean_dif_J_sc=mean(abs(f1-J_sc)./J_sc);
NRMSE_J_sc= sqrt(sum((f1-J_sc).^2)/(N_curves*mean(J_sc)));
figuretext=sprintf('Fit of Jsc\nNRMSE=%.2f\nRMD=%.4f',NRMSE_J_sc,Rel_mean_dif_J_sc);
text(300, 0.15*f(1), figuretext)

%% Fit Voc
V_oc=par_T(:,2)';
po_V_oc=polyfit(T_levels,V_oc,3);

%Plot
figure;
f=polyval(po_V_oc,T_1K);
hold on
plot(T_1K,f,T_levels,V_oc,'o');
title('Open Circuit Voltage vs. Irradiance');
box on
xlabel('Irradiance (W/m^2)');
ylabel('V_O_C (V)');
hold off

%Determine goodnes of fit
f1=(polyval(po_V_oc,T_levels));
Rel_mean_dif_V_oc=mean(abs(f1-V_oc)./V_oc);
NRMSE_V_oc= sqrt(sum((f1-V_oc).^2)/(N_curves*mean(V_oc)));
figuretext=sprintf('Fit of Voc\nNRMSE=%.2f\nRMD=%.4f',NRMSE_V_oc,Rel_mean_dif_V_oc);
text(300, 1.05*f(1), figuretext)




%% Fit Pmmp
Jmpp=par_T(:,8)'*(246/10000); %Convert to A
Vmpp=par_T(:,9)';

P_mmp=Jmpp.*Vmpp;
po_P_mmp=polyfit(T_levels,P_mmp,3);

%Plot
figure;
f=polyval(po_P_mmp,T_1K);
hold on
plot(T_1K,f,T_levels,P_mmp,'o');
title('Maximum power point  vs. Irradiance');
box on
xlabel('Temperature (W/m^2)');
ylabel('P_{mpp} (W)');
hold off

%Determine goodnes of fit
f1=(polyval(po_P_mmp,T_levels));
Rel_mean_dif_P_mmp=mean(abs(f1-P_mmp)./P_mmp);
NRMSE_P_mmp= sqrt(sum((f1-P_mmp).^2)/(N_curves*mean(P_mmp)));
figuretext=sprintf('Fit of Voc\nNRMSE=%.2f\nRMD=%.4f',NRMSE_P_mmp,Rel_mean_dif_P_mmp);
text(300, 1.05*f(1), figuretext)


end