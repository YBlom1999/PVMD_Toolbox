clear;
%Define irradiances of the IV curves
G=[50,100,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500];
numOfVariations=length(G);
T=298.15;

%Line for starting and ending the import
startLine=7;
endLine=157;
numOfLines=endLine-startLine+1;

Voltage=zeros(numOfLines,numOfVariations);
Current=zeros(numOfLines,numOfVariations);
par_G=zeros(numOfVariations,9);

%Set file names
for i=1:numOfVariations
    filename=join(['JV_pk_HZB_2T_Top_v2_25deg_',num2str(G(i)),'Wm2.dat']);
   
    [Voltage(:,i), Current(:,i)] = importIVCurve(filename, [startLine, endLine]);

end

%Convert to A/m2
%Current=10*Current;
Current=Current*-1.0;


for j=1:numOfVariations
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
    f_rp=fit(V(ind_Jsc:ind_Jsc+75),J(ind_Jsc:ind_Jsc+75),'poly1');
    coeff_f_rp=coeffvalues(f_rp);
    Rp=-1.0/coeff_f_rp(1);
    
    
    rp=Rp/(Voc/Jsc);
    Vth=(1.3806*10^-23)*T/(1.602*10^-19);
    % x(1)=voc - x(2)=rs
    FF0=@(x) ((x(1) - log(x(1)+0.72))/(x(1)+1));
    FFs=@(x) FF0(x)*(1-1.1*x(2))+(x(2)^2)/5.4;
    FFsh=@(x) FFs(x)*(1-((x(1)+0.7)/x(1))*(FFs(x)/rp));
    fill=@(x) abs(FF - FF0(x));
    
    
    if j==1
        X0=[Jsc/Voc];
    else
        X0=[voc];
    end
    [voc,err]=fminsearch(fill,X0);
    %
    fill=@(x) abs(FF - FFsh(x));
    X0=[voc 1e-4];
    [voc,err]=fminsearch(fill,X0);
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
    [J0,err]=fminsearch(Err,X0,options);
    
    
    par_G(j,:)=[Jsc, Voc, Rp, Rs, J0, n, FF, Jmpp, Vmpp];
    
    %%PLot
    clf;
    hold on
    Jmd=Jsc - J0.*(exp((V+J*Rs)./(Vth*n))-1)-(V+Rs.*J)./Rp;
    plot(V(1:k),J(1:k),V(1:k),Jmd(1:k))
    ylabel('Current density [A/m2]')
    xlabel('Voltage [V]')
    box on
    
    Rel_mean_dif=mean(abs(Jmd(1:k)-J(1:k))./J(1:k));
    NRMSE= sqrt(sum((Jmd(1:k)-J(1:k)).^2)/(k*mean(J(1:k))));
    figuretext=sprintf('JV Fit G=%d\nNRMSE=%.2f\nRMD=%.4f',G(j),NRMSE,Rel_mean_dif);
    text(0.1, 0.3*J(1), figuretext)
    hold off
    
    figurename=join(['JV_Fit_G',num2str(G(j)),'.fig']);
    savefig(figurename);
    
end


%% Fit Jsc
J=par_G(:,1)';
po_Jph_G=polyfit(G,J,1);
G_1_W_m2=min(G):max(G);


%plot
clf;
J_values=polyval(po_Jph_G,G_1_W_m2);
hold on
plot(G_1_W_m2,J_values/10,'LineWidth',2);
plot(G,J/10,'o','MarkerSize',8);
ax = gca;
ax.FontSize = 13; 
xlim([min(G) max(G)]);
hold off
box on
xlabel('Irradiance (W/m^2)');
ylabel('Photo-generated Current Density [mA/cm^2]');

%Determine goodnes of fit
f1=polyval(po_Jph_G,G);
Rel_mean_dif_J=mean(abs(f1-J)./J);
NRMSE_J= sqrt(sum((f1-J).^2)/(numOfVariations*mean(J)));
display(Rel_mean_dif_J);
display(NRMSE_J);
figuretext=sprintf('Fit of Jsc\nNRMSE=%.2f\nRMD=%.4f',NRMSE_J,Rel_mean_dif_J);
text(200, 0.15*J_values(1), figuretext)
savefig('Fit_J_G.fig');

%the Jsc which corresponds to G at every level, will be used in polyfit, 
%since in SimulatedIV j_ph is know but not G




%% Fit shunt resistance
shunt=par_G(:,3)';
R_sh_J=polyfit(J,log(shunt),4);

%plot
clf;

f=polyval(R_sh_J,J_values);
f=exp(f);
hold on
plot(G_1_W_m2,f/1e2,'LineWidth',2)
plot(G,shunt/1e2,'o','MarkerSize',8)
% grid on
%title('Shunt Resistance vs. Irradiance')
ylabel('Shunt Resistance [MOhomcm^2]')
xlabel('Irradiance Intensity [W/m^2]')
hold off
box on
ax = gca;
ax.FontSize = 13; 
%ylim([0 6])
xlim([min(G) max(G)])
% RMSE
f1=exp(polyval(R_sh_J,J));
Rel_mean_dif_Rsh=mean(abs(f1-shunt)./shunt);
NRMSE_Rsh= sqrt(sum((f1-shunt).^2)/(numOfVariations*mean(shunt)));
display(Rel_mean_dif_Rsh);
display(NRMSE_Rsh);
figuretext=sprintf('Fit of Rsh\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Rsh,Rel_mean_dif_Rsh);
text(200, 0.8*f(1)/1e2, figuretext)

savefig('Fit_Rsh_J.fig');


%% Fit series resitance
series=par_G(:,4)';
R_s_J=polyfit(J,log(series),3);

%plot
clf;
f=polyval(R_s_J,J_values);
f=exp(f);
hold on
plot(G_1_W_m2,f*1e4,'LineWidth',2);
plot(G,series*1e4,'o','MarkerSize',8);
ylabel('Series Resistance [Ohm.cm^2]');
xlabel('Irradiance Intensity [W/m^2]');
ax = gca;
ax.FontSize = 13; 
xlim([min(G) max(G)])
%ylim([0 0.017])
box on
hold off

% RMSE
f1=exp(polyval(R_s_J,J));
Rel_mean_dif_Rs=mean(abs(f1-series)./series);
NRMSE_Rs= sqrt(sum((f1-series).^2)/(numOfVariations*mean(series)));
display(Rel_mean_dif_Rs);
display(NRMSE_Rs);
figuretext=sprintf('Fit of Rs\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Rs,Rel_mean_dif_Rs);
text(200, 0.8*f(1)*1e4, figuretext)
savefig('Fit_Rs_J.fig');

%% Fit J0
sat=par_G(:,5)';
J0_J=polyfit(J,log(sat),3);

%Plot
clf;
f=exp(polyval(J0_J,J_values));
hold on
plot(G_1_W_m2,f*1e8,'LineWidth',2);
plot(G,sat*1e8,'o','MarkerSize',8);
xlabel('Irradiance Intensity [W/m^2]');
ylabel('Saturation Current Density [pA/cm^2]');
ax = gca;
ax.FontSize = 13; 
ax.YScale='log';
xlim([min(G) max(G)])
box on
hold off

% RMSE
f1=exp(polyval(J0_J,J));
Rel_mean_dif_J0=mean(abs(f1-sat)./sat);
NRMSE_J0= sqrt(sum((f1-sat).^2)/(numOfVariations*mean(sat)));
display(Rel_mean_dif_J0);
display(NRMSE_J0);
figuretext=sprintf('Fit of J0\nNRMSE=%.2f\nRMD=%.4f',NRMSE_J0,Rel_mean_dif_J0);
text(200, 0.8*f(1)*1e8, figuretext)
savefig('Fit_J0_J.fig');


%% Fit idiality n
n=par_G(:,6)';
n_J=polyfit(J,n,3);

% Plot
clf;
f=polyval(n_J,J_values);
hold on
plot(G_1_W_m2,f,'LineWidth',2);
plot(G,n,'o','MarkerSize',8);
hold off
ax = gca;
ax.FontSize = 13; 
%ylim([0.8 1.2]);
xlim([min(G) max(G)]);
box on
xlabel('Irradiance Intensity [W/m^2]');
ylabel('Ideality Factor');

%RMSE
f1=polyval(n_J,J);
Rel_mean_dif_n=mean(abs(f1-n)./n);
NRMSE_n= sqrt(sum((f1-n).^2)/(numOfVariations*mean(n)));
display(Rel_mean_dif_n);
display(NRMSE_n);
figuretext=sprintf('Fit of idialty factor\nNRMSE=%.2f\nRMD=%.4f',NRMSE_n,Rel_mean_dif_n);
text(800, f(1), figuretext)
savefig('Fit_n_J.fig');

%% Fit fill factor
Fill=par_G(:,7)';
po=polyfit(J,Fill,3);

%plot
clf;
hold on
f=polyval(po,J_values);
plot(G_1_W_m2,f,G,Fill,'o');
hold off

box on
xlabel('Irradiance (W/m^2)');
ylabel('FF');

%RMSE
f1=polyval(po,J);
Rel_mean_dif_Fill=mean(abs(f1-Fill)./Fill);
NRMSE_Fill= sqrt(sum((f1-Fill).^2)/(numOfVariations*mean(Fill)));
display(Rel_mean_dif_Fill);
display(NRMSE_Fill);
figuretext=sprintf('Fit of fill factor\nNRMSE=%.2f\nRMD=%.4f',NRMSE_Fill,Rel_mean_dif_Fill);
text(200, 1.01*f(1), figuretext);
savefig('Fit_Fill_J.fig');

%% Fit Voc 
V_oc=par_G(:,2)';
po_V_oc=polyfit(J,V_oc,3);

%Plot
clf;
f=polyval(po_V_oc,J_values);
hold on
plot(G_1_W_m2,f,G,V_oc,'o');
title('Open Circuit Voltage vs. Irradiance');
box on
xlabel('Irradiance (W/m^2)');
ylabel('V_O_C (V)');
hold off

%Determine goodnes of fit
f1=(polyval(po_V_oc,J));
Rel_mean_dif_V_oc=mean(abs(f1-V_oc)./V_oc);
NRMSE_V_oc= sqrt(sum((f1-V_oc).^2)/(numOfVariations*mean(V_oc)));
display(Rel_mean_dif_V_oc);
display(NRMSE_V_oc);
figuretext=sprintf('Fit of Voc\nNRMSE=%.2f\nRMD=%.4f',NRMSE_V_oc,Rel_mean_dif_V_oc);
text(200, 1.05*f(1), figuretext)
savefig('Fit_V_oc_J.fig');

save 1Dode_J.mat J0_J R_s_J R_sh_J n_J;

