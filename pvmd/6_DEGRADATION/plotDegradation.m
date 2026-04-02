function plotDegradation(CELL_output,CELL_output_new,ELECTRIC_output,ELECTRIC_output_new,Time,k_dis,k_mois,k_LID,k_TC,k_total,numCells,Rcon_init,Rcon_new,CONSTANTS)

k_b = CONSTANTS.k_b;
q = CONSTANTS.q;
T_STC = CONSTANTS.T_STC;
Vth=k_b*T_STC/q;

%Plot Degradation curves
k_dis_avg = mean(reshape(k_dis,[24*5,length(k_dis)/24/5]));
k_mois_avg = mean(reshape(k_mois,[24*5,length(k_mois)/24/5]));
k_LID_avg = mean(reshape(k_LID,[24*5,length(k_LID)/24/5]));
k_TC_avg = mean(reshape(k_TC,[24*5,length(k_TC)/24/5]));
k_total_avg = mean(reshape(k_total,[24*5,length(k_total)/24/5]));
Time_avg = 24*5:24*5:Time(end);

figure
hold on; box on;
plot(Time_avg/8760,k_dis_avg)
plot(Time_avg/8760,k_mois_avg)
plot(Time_avg/8760,k_LID_avg)
plot(Time_avg/8760,k_TC_avg)
plot(Time_avg/8760,k_total_avg,'color','k')
xlabel('Time [year]')
ylabel('Degradation rate [-]')
legend('Discoloration','Moisture induced','LID','Thermal cycling')
ylim([0,1.2*max(k_total_avg)]);
ax = gca;
ax.FontSize = 15;

%plot EQE and IV curve
Type = CELL_output.TYPE;
wav = CELL_output.CELL_FRONT.wav;
Voltage = 0:0.01:1.4;

if strcmp(Type,'Tan')
    EQE1_init = CELL_output.CELL_FRONT.RAT(:,1,2);
    EQE1_new = CELL_output_new.CELL_FRONT.RAT(:,1,2);

    EQE2_init = CELL_output.CELL_FRONT.RAT(:,1,3);
    EQE2_new = CELL_output_new.CELL_FRONT.RAT(:,1,3);

    Iph1_init = ELECTRIC_output.Parameters_STC_1(1,1);
    Rs1_init = ELECTRIC_output.Parameters_STC_1(1,2);
    Rsh1_init = ELECTRIC_output.Parameters_STC_1(1,3);
    n1_init = ELECTRIC_output.Parameters_STC_1(1,4);
    I01_init = ELECTRIC_output.Parameters_STC_1(1,5);
   
    Iph2_init = ELECTRIC_output.Parameters_STC_2(1,1);
    Rs2_init = ELECTRIC_output.Parameters_STC_2(1,2);
    Rsh2_init = ELECTRIC_output.Parameters_STC_2(1,3);
    n2_init = ELECTRIC_output.Parameters_STC_2(1,4);
    I02_init = ELECTRIC_output.Parameters_STC_2(1,5);

    Current_axes = 0:0.01:max([Iph1_init,Iph2_init])*1.2;

    z=(Rs1_init*I01_init/(n1_init*Vth*(1+Rs1_init/Rsh1_init)))*exp((Rs1_init*(Iph1_init+I01_init)+Voltage)./(n1_init*Vth*(1+Rs1_init/Rsh1_init)));
    Current=(Iph1_init+I01_init-Voltage/(Rsh1_init))/(1+Rs1_init/Rsh1_init)-lambertw(z).*(n1_init*Vth)/Rs1_init;
    Voltage1_init = interp1(Current,Voltage,Current_axes,'linear','extrap');

    z=(Rs2_init*I02_init/(n2_init*Vth*(1+Rs2_init/Rsh2_init)))*exp((Rs2_init*(Iph2_init+I02_init)+Voltage)./(n2_init*Vth*(1+Rs2_init/Rsh2_init)));
    Current=(Iph2_init+I02_init-Voltage/(Rsh2_init))/(1+Rs2_init/Rsh2_init)-lambertw(z).*(n2_init*Vth)/Rs2_init;
    Voltage2_init = interp1(Current,Voltage,Current_axes,'linear','extrap');

    Voltage_init = Voltage1_init+Voltage2_init;
    Voltage_init = (Voltage_init-Current_axes*Rcon_init)*numCells;

    Iph1_new = ELECTRIC_output_new.Parameters_STC_1(1,1);
    Rs1_new = ELECTRIC_output_new.Parameters_STC_1(1,2);
    Rsh1_new = ELECTRIC_output_new.Parameters_STC_1(1,3);
    n1_new = ELECTRIC_output_new.Parameters_STC_1(1,4);
    I01_new = ELECTRIC_output_new.Parameters_STC_1(1,5);
   
    Iph2_new = ELECTRIC_output_new.Parameters_STC_2(1,1);
    Rs2_new = ELECTRIC_output_new.Parameters_STC_2(1,2);
    Rsh2_new = ELECTRIC_output_new.Parameters_STC_2(1,3);
    n2_new = ELECTRIC_output_new.Parameters_STC_2(1,4);
    I02_new = ELECTRIC_output_new.Parameters_STC_2(1,5);

    z=(Rs1_new*I01_new/(n1_new*Vth*(1+Rs1_new/Rsh1_new)))*exp((Rs1_new*(Iph1_new+I01_new)+Voltage)./(n1_new*Vth*(1+Rs1_new/Rsh1_new)));
    Current=(Iph1_new+I01_new-Voltage/(Rsh1_new))/(1+Rs1_new/Rsh1_new)-lambertw(z).*(n1_new*Vth)/Rs1_new;
    Voltage1_new = interp1(Current,Voltage,Current_axes,'linear','extrap');

    z=(Rs2_new*I02_new/(n2_new*Vth*(1+Rs2_new/Rsh2_new)))*exp((Rs2_new*(Iph2_new+I02_new)+Voltage)./(n2_new*Vth*(1+Rs2_new/Rsh2_new)));
    Current=(Iph2_new+I02_new-Voltage/(Rsh2_new))/(1+Rs2_new/Rsh2_new)-lambertw(z).*(n2_new*Vth)/Rs2_new;
    Voltage2_new = interp1(Current,Voltage,Current_axes,'linear','extrap');

    Voltage_new = Voltage1_new+Voltage2_new;
    Voltage_new = (Voltage_new-Current_axes*Rcon_new)*numCells;

    fig = figure;
    fig.Position = fig.Position + [-300,0,600,0];
    subplot(1,2,1)
    hold on; box on;
    plot(wav*1e3,EQE1_init,'LineWidth',2,'color','b')
    plot(wav*1e3,EQE1_new,'--','LineWidth',2,'color','b');
    plot(wav*1e3,EQE2_init,'LineWidth',2,'color','r')
    plot(wav*1e3,EQE2_new,'--','LineWidth',2,'color','r');
    xlabel('Wavelength [nm]')
    ylabel('EQE [-]')
    xlim([wav(1)*1e3,wav(end)*1e3])
    ylim([0, 1])
    legend('Top initial','Top after','Bottom initial','Bottom after','Location','south')
    ax = gca;
    ax.FontSize = 15;

    subplot(1,2,2)
    hold on;
    box on;
    plot(Voltage_init,Current_axes,'LineWidth',2)
    plot(Voltage_new,Current_axes,'LineWidth',2)
    xlabel('Voltage [V]')
    ylabel('Current [A]')
    xlim([0,1.2*max(Voltage_init)])
    legend('Initial','After','Location','west')
    ax = gca;
    ax.FontSize = 15;

    figure
    hold on;
    box on;
    plot(Voltage1_init,Current_axes,'LineWidth',2,'color','b')
    plot(Voltage1_new,Current_axes,'--','LineWidth',2,'color','b');
    plot(Voltage2_init,Current_axes,'LineWidth',2,'color','r')
    plot(Voltage2_new,Current_axes,'--','LineWidth',2,'color','r');
    xlabel('Voltage [V]')
    ylabel('Current [A]')
    xlim([0,1.5])
    legend('Top initial','Top after','Bottom initial','Bottom after','Location','southwest')
    ax = gca;
    ax.FontSize = 15;

else

    EQE_init = CELL_output.CELL_FRONT.RAT(:,1,2);
    EQE_new = CELL_output_new.CELL_FRONT.RAT(:,1,2);

    Iph_init = ELECTRIC_output.Parameters_STC_1(1,1);
    Rs_init = ELECTRIC_output.Parameters_STC_1(1,2);
    Rsh_init = ELECTRIC_output.Parameters_STC_1(1,3);
    n_init = ELECTRIC_output.Parameters_STC_1(1,4);
    I0_init = ELECTRIC_output.Parameters_STC_1(1,5);

    Current_axes = 0:0.01:Iph_init*1.2;

    z=(Rs_init*I0_init/(n_init*Vth*(1+Rs_init/Rsh_init)))*exp((Rs_init*(Iph_init+I0_init)+Voltage)./(n_init*Vth*(1+Rs_init/Rsh_init)));
    Current=(Iph_init+I0_init-Voltage/(Rsh_init))/(1+Rs_init/Rsh_init)-lambertw(z).*(n_init*Vth)/Rs_init;
    Voltage_init = interp1(Current,Voltage,Current_axes,'linear','extrap');
    Voltage_init = (Voltage_init-Current_axes*Rcon_init)*numCells;

    Iph_new = ELECTRIC_output_new.Parameters_STC_1(1,1);
    Rs_new = ELECTRIC_output_new.Parameters_STC_1(1,2);
    Rsh_new = ELECTRIC_output_new.Parameters_STC_1(1,3);
    n_new = ELECTRIC_output_new.Parameters_STC_1(1,4);
    I0_new = ELECTRIC_output_new.Parameters_STC_1(1,5);

    z=(Rs_new*I0_new/(n_new*Vth*(1+Rs_new/Rsh_new)))*exp((Rs_new*(Iph_new+I0_new)+Voltage)./(n_new*Vth*(1+Rs_new/Rsh_new)));
    Current=(Iph_new+I0_new-Voltage/(Rsh_new))/(1+Rs_new/Rsh_new)-lambertw(z).*(n_new*Vth)/Rs_new;
    Voltage_new = interp1(Current,Voltage,Current_axes,'linear','extrap');
    Voltage_new = (Voltage_new-Current_axes*Rcon_new)*numCells;

    fig = figure;
    fig.Position = fig.Position + [-300,0,600,0];
    subplot(1,2,1)
    hold on; box on;
    plot(wav*1e3,EQE_init,'LineWidth',2)
    plot(wav*1e3,EQE_new,'LineWidth',2);
    xlabel('Wavelength [nm]')
    ylabel('EQE [-]')
    xlim([wav(1)*1e3,wav(end)*1e3])
    ylim([0, 1])
    legend('Initial','After','Location','south')
    ax = gca;
    ax.FontSize = 15;

    subplot(1,2,2)
    hold on;
    box on;
    plot(Voltage_init,Current_axes,'LineWidth',2)
    plot(Voltage_new,Current_axes,'LineWidth',2)
    xlabel('Voltage [V]')
    ylabel('Current [A]')
    xlim([0,1.2*max(Voltage_init)])
    legend('Initial','After','Location','west')
    ax = gca;
    ax.FontSize = 15;


end

end