function [J_sc_J,J0_J,R_s_J,R_sh_J,n_J] = Fit_Irr_curves(Voltage, Current,Jph,Irr_levels,T_fixed)
%Fit_Irr_curves Calculates the dependencies of the equivalent circuit
%parameters w.r.t. the absorbed current
%
% This function reads the IV curves, fits it with the one-diode equivalent
% circuit and makes a fit w.r.t. the absorbed current
%
% Parameters
% ----------
% Voltage: double
%   The voltages of all the IV curves
% Current: double
%   The currents of all the IV curves
% Jph: double
%   The absorbed photogenerated current of all the IV curves
% Irr_levels : double
%   The irradiance levels for which IV curves are measured/simulated.
% T_fixed: double
%   The temperature for which all measurements are done.
%
% Returns
% -------
% J_sc_J : double
%   The values of the fit for the short circuit current
% J_0_J : double
%   The values of the fit for the saturation current
% R_s_J : double
%   The values of the fit for the series resistance
% R_sh_J : double
%   The values of the fit for the shunt resistance
% n_J : double
%   The values of the fit for the ideality factor
%
% Developed by Malte Vogt and Youri Blom

N_curves = size(Voltage,2);
T_fixed=T_fixed + 273.15;
par_G=zeros(N_curves,9);

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
    %f_rp=fit(V(ind_Jsc-5:ind_Jsc+5),J(ind_Jsc-5:ind_Jsc+5),'poly1');
    f_rp=fit(V(ind_Jsc:ind_Jsc+5),J(ind_Jsc:ind_Jsc+5),'poly1');
    coeff_f_rp=coeffvalues(f_rp);
    Rp=-1.0/coeff_f_rp(1);
    
    
    rp=Rp/(Voc/Jsc);
    Vth=(1.3806*10^-23)*T_fixed/(1.602*10^-19);
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
    
    
    par_G(j,:)=[Jsc, Voc, Rp, Rs, J0, n, FF, Jmpp, Vmpp];
    
    %%PLot
    Jmd=Jsc - J0.*(exp((V+J*Rs)./(Vth*n))-1)-(V+Rs.*J)./Rp;
    plot(V,J,'color',color(j,:))
    plot(V,Jmd,'--','color',color(j,:))

       
    
end
xlabel('Voltage [V]')
ylabel('Current density [A/m2]')
xlim([0, max(max(Voltage))])
ylim([0, 1.05*max(Current(1,:))]);

%% Fit Jph
Jsc = par_G(:,1)';
po_Jsc_G=polyfit(Irr_levels,Jph,1);
po_Jph_G=polyfit(Irr_levels,Jph,1);
G_1_W_m2=min(Irr_levels):max(Irr_levels);


%plot
figure;
Jsc_values=polyval(po_Jsc_G,G_1_W_m2);
Jph_values=polyval(po_Jph_G,G_1_W_m2);
hold on
plot(G_1_W_m2,Jph_values/10,'LineWidth',2);
plot(Irr_levels,Jph/10,'o','MarkerSize',8);
ax = gca;
ax.FontSize = 13; 
xlim([min(Irr_levels) max(Irr_levels)]);
hold off
box on
xlabel('Irradiance (W/m^2)');
ylabel('Photo-generated Current Density [mA/cm^2]');






%% Fit shunt resistance
shunt=par_G(:,3)';
R_sh_J=polyfit(Jph,log(shunt),4);

%plot
figure;

f=polyval(R_sh_J,Jph_values);
f=exp(f);
hold on
plot(G_1_W_m2,f/1e2,'LineWidth',2)
plot(Irr_levels,shunt/1e2,'o','MarkerSize',8)
% grid on
%title('Shunt Resistance vs. Irradiance')
ylabel('Shunt Resistance [MOhomcm^2]')
xlabel('Irradiance Intensity [W/m^2]')
hold off
box on
ax = gca;
ax.FontSize = 13; 
%ylim([0 6])
xlim([min(Irr_levels) max(Irr_levels)])
% RMSE
f1=exp(polyval(R_sh_J,Jph));
Rel_mean_dif_Rsh=mean(abs(f1-shunt)./shunt);
NRMSE_Rsh= sqrt(sum((f1-shunt).^2)/(N_curves*mean(shunt)));
%display(Rel_mean_dif_Rsh);
%display(NRMSE_Rsh);
figuretext=sprintf('Fit of Rsh\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Rsh,Rel_mean_dif_Rsh);
text(200, 0.8*f(1)/1e2, figuretext)



%% Fit series resitance
series=par_G(:,4)';
R_s_J=polyfit(Jph,log(series),3);

%plot
figure;
f=polyval(R_s_J,Jph_values);
f=exp(f);
hold on
plot(G_1_W_m2,f*1e4,'LineWidth',2);
plot(Irr_levels,series*1e4,'o','MarkerSize',8);
ylabel('Series Resistance [Ohm.cm^2]');
xlabel('Irradiance Intensity [W/m^2]');
ax = gca;
ax.FontSize = 13; 
xlim([min(Irr_levels) max(Irr_levels)])
%ylim([0 0.017])
box on
hold off

% RMSE
f1=exp(polyval(R_s_J,Jph));
Rel_mean_dif_Rs=mean(abs(f1-series)./series);
NRMSE_Rs= sqrt(sum((f1-series).^2)/(N_curves*mean(series)));
%display(Rel_mean_dif_Rs);
%display(NRMSE_Rs);
figuretext=sprintf('Fit of Rs\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Rs,Rel_mean_dif_Rs);
text(200, 0.8*f(1)*1e4, figuretext)

%% Fit J0
sat=par_G(:,5)';
J0_J=polyfit(Jph,log(sat),3);

%Plot
figure;
f=exp(polyval(J0_J,Jsc_values));
hold on
plot(G_1_W_m2,f*1e8,'LineWidth',2);
plot(Irr_levels,sat*1e8,'o','MarkerSize',8);
xlabel('Irradiance Intensity [W/m^2]');
ylabel('Saturation Current Density [pA/cm^2]');
ax = gca;
ax.FontSize = 13; 
ax.YScale='log';
xlim([min(Irr_levels) max(Irr_levels)])
box on
hold off

% RMSE
f1=exp(polyval(J0_J,Jsc));
Rel_mean_dif_J0=mean(abs(f1-sat)./sat);
NRMSE_J0= sqrt(sum((f1-sat).^2)/(N_curves*mean(sat)));
%display(Rel_mean_dif_J0);
%display(NRMSE_J0);
figuretext=sprintf('Fit of J0\nNRMSE=%.2f\nRMD=%.4f',NRMSE_J0,Rel_mean_dif_J0);
text(200, 0.8*f(1)*1e8, figuretext)

%% Fit idiality n
n=par_G(:,6)';
n_J=polyfit(Jph,n,4);

% Plot
figure;
f=polyval(n_J,Jph_values);
hold on
plot(G_1_W_m2,f,'LineWidth',2);
plot(Irr_levels,n,'o','MarkerSize',8);
hold off
ax = gca;
ax.FontSize = 13;
%ylim([0.8 1.2]);
xlim([min(Irr_levels) max(Irr_levels)]);
box on
xlabel('Irradiance Intensity [W/m^2]');
ylabel('Ideality Factor');

%RMSE
f1=polyval(n_J,Jph);
Rel_mean_dif_n=mean(abs(f1-n)./n);
NRMSE_n= sqrt(sum((f1-n).^2)/(N_curves*mean(n)));
%display(Rel_mean_dif_n);
%display(NRMSE_n);
figuretext=sprintf('Fit of idialty factor\nNRMSE=%.2f\nRMD=%.4f',NRMSE_n,Rel_mean_dif_n);
text(300, f(1), figuretext)

%% Fit Short circuit current
J_sc_J=polyfit(Jph,Jsc,3);

% Plot
figure;
f=polyval(J_sc_J,Jph_values);
hold on
plot(G_1_W_m2,f,'LineWidth',2);
plot(Irr_levels,Jsc,'o','MarkerSize',8);
hold off
ax = gca;
ax.FontSize = 13; 
%ylim([0.8 1.2]);
xlim([min(Irr_levels) max(Irr_levels)]);
box on
xlabel('Irradiance Intensity [W/m^2]');
ylabel('J_{sc} [A/m^2]');

%RMSE
f1=polyval(J_sc_J,Jph);
Rel_mean_dif_Jsc=mean(abs(f1-n)./n);
NRMSE_Jsc= sqrt(sum((f1-n).^2)/(N_curves*mean(n)));
%display(Rel_mean_dif_Jsc);
%display(NRMSE_Jsc);
figuretext=sprintf('Fit of idialty factor\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Jsc,Rel_mean_dif_Jsc);
text(800, f(1), figuretext)

%% Fit fill factor
Fill=par_G(:,7)';
po=polyfit(Jph,Fill,3);

%plot
figure
hold on
f=polyval(po,Jph_values);
plot(G_1_W_m2,f,Irr_levels,Fill,'o');
hold off

box on
xlabel('Irradiance (W/m^2)');
ylabel('FF');

%RMSE
f1=polyval(po,Jph);
Rel_mean_dif_Fill=mean(abs(f1-Fill)./Fill);
NRMSE_Fill= sqrt(sum((f1-Fill).^2)/(N_curves*mean(Fill)));
%display(Rel_mean_dif_Fill);
%display(NRMSE_Fill);
figuretext=sprintf('Fit of fill factor\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Fill,Rel_mean_dif_Fill);
text(200, 1.01*f(1), figuretext);

%% Fit Voc 
V_oc=par_G(:,2)';
po_V_oc=polyfit(Jph,V_oc,3);

%Plot
figure
f=polyval(po_V_oc,Jph_values);
hold on
plot(G_1_W_m2,f,Irr_levels,V_oc,'o');
title('Open Circuit Voltage vs. Irradiance');
box on
xlabel('Irradiance (W/m^2)');
ylabel('V_O_C (V)');
hold off

%Determine goodnes of fit
f1=(polyval(po_V_oc,Jph));
Rel_mean_dif_V_oc=mean(abs(f1-V_oc)./V_oc);
NRMSE_V_oc= sqrt(sum((f1-V_oc).^2)/(N_curves*mean(V_oc)));
%display(Rel_mean_dif_V_oc);
%display(NRMSE_V_oc);
figuretext=sprintf('Fit of Voc\nNRMSE=%.2f\nRMD=%.4f',NRMSE_V_oc,Rel_mean_dif_V_oc);
text(200, 1.05*f(1), figuretext)


end