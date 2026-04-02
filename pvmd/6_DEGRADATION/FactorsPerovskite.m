function [factors] = FactorsPerovskite(Parameters,k_needed,Scenario,T)
%Degradation_single file for the degradation module of tandem modules
%
% This function calculates the degradation rate of a tandem module
% 
%
% Parameters
% ----------
% Parameters : double
%   Parameters of the equivalent circuit
% k_needed : double
%   The degradation rate that is needed.
% Scenario : double
%   The degradation scenario
%
% Returns
% -------
% factors : double
%   The factors by which the parameters should be updated
%
% Developed by by Youri Blom

K=1.3806e-23;
q=1.6022e-19;
Vth=K*T/q;
V = 0:0.01:1.3;

Iph = Parameters(1);
Rs = Parameters(2);
Rsh = Parameters(3);
n = Parameters(4);
I0 = Parameters(5);

z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+V)./(n*Vth*(1+Rs/Rsh)));
J=(Iph+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;

P_orig = max(J.*V);
if P_orig == 0
    factors = [1,1,1,1,1];
    return
end
P_test = zeros(500,1);
if Scenario == 1
    for i = 1:500
        Iph_test = Iph*(1-(i-1)/500);
        z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph_test+I0)+V)./(n*Vth*(1+Rs/Rsh)));
        J=(Iph_test+I0-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;
        P_test(i) = max(J.*V);
        if 1-P_test(i)/P_orig > k_needed; break; end
    
    end
    P_loss = 1-P_test./P_orig;
    [~,ind] = unique(P_loss);
    i_needed = interp1(P_loss(ind),ind,k_needed,'linear','extrap');
    Iph_new = Iph*(1-(i_needed-1)/500);
    n_new = n;
    Rsh_new = Rsh;
    Rs_new = Rs;
    I0_new = I0;
elseif Scenario == 2
    for i = 1:500
        I0_test = exp(log(I0)+(i-1)/15);
        z=(Rs*I0_test/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0_test)+V)./(n*Vth*(1+Rs/Rsh)));
        J=(Iph+I0_test-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;
        P_test(i) = max(J.*V);
        if 1-P_test(i)/P_orig > k_needed; break; end
    end
    P_loss = 1-P_test./P_orig;
    
    [~,ind] = unique(P_loss);
    i_needed = interp1(P_loss(ind),ind,k_needed,'linear','extrap');
    Iph_new = Iph;
    n_new = n;
    Rsh_new = Rsh;
    Rs_new = Rs;
    I0_new = exp(log(I0)+(i_needed-1)/15);
    
elseif Scenario == 3
    Voc_orig = n*Vth*log(Iph/I0);
    for i = 1:500
        n_test = n*(1+(i-1)/25);
        I0_test = Iph/(exp(Voc_orig/n_test/Vth));
        z=(Rs*I0_test/(n_test*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0_test)+V)./(n_test*Vth*(1+Rs/Rsh)));
        J=(Iph+I0_test-V/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n_test*Vth)/Rs;
        P_test(i) = max(J.*V);
        if 1-P_test(i)/P_orig > k_needed; break; end
        
    end
    P_loss = 1-P_test./P_orig;
    
    [~,ind] = unique(P_loss);
    i_needed = interp1(P_loss(ind),ind,k_needed,'linear','extrap');
    Iph_new = Iph;
    n_new = n*(1+(i_needed-1)/25);
    Rsh_new = Rsh;
    Rs_new = Rs;
    I0_new = Iph/(exp(Voc_orig/n_new/Vth));

elseif Scenario == 4
    Rs_final = 1;
    Rsh_final = 1;
    for i = 1:500
        step = (i-1)/length(P_test);
        Rs_test = exp((1-step)*log(Rs)+step*log(Rs_final));%1./(1/Rs*(1-step)+1/Rs_final*step);
        Rsh_test = 1/(1/Rsh*(1-step)+1/Rsh_final*step);%Rsh*(1-step) + Rsh_final*step;
        z=(Rs_test*I0/(n*Vth*(1+Rs_test/Rsh_test)))*exp((Rs_test*(Iph+I0)+V)./(n*Vth*(1+Rs_test/Rsh_test)));
        J=(Iph+I0-V/(Rsh_test))/(1+Rs_test/Rsh_test)-lambertw(z).*(n*Vth)/Rs_test;
        P_test(i) = max(J.*V);
        if 1-P_test(i)/P_orig > k_needed; break; end
    end
    P_loss = 1-P_test./P_orig;
    
    [~,ind] = unique(P_loss);
    i_needed = interp1(P_loss(ind),ind,k_needed,'linear','extrap');
    step_needed = (i_needed-1)/length(P_test);
    Iph_new = Iph;
    n_new = n;
    Rsh_new = 1/(1/Rsh*(1-step_needed)+1/Rsh_final*step_needed);%Rsh*(1-step_needed) + Rsh_final*step_needed;
    Rs_new = exp((1-step_needed)*log(Rs)+step_needed*log(Rs_final));%1./(1/Rs*(1-step_needed)+1/Rs_final*step_needed);
    I0_new = I0;
end

factor_Iph = Iph_new/Iph;
factor_Rs = Rs_new/Rs;
factor_Rsh = Rsh_new/Rsh;
factor_n = n_new/n;
factor_I0 = log10(I0_new)/log10(I0);

factors = [factor_Iph,factor_Rs,factor_Rsh,factor_n,factor_I0];

end