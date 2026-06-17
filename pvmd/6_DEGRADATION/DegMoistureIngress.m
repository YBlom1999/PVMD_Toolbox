function [k_pred] = DegMoistureIngress(TOOLBOX_input,T,Ea,A,C,Time)
%DegMoistureIngress Calculates the degradation due to moisture ingress
%
% This function calculates the degradation rate caused by the moisture
% ingress into the PV module
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% T : double
%   Temperature of the module
% Ea : double
%   The activation energy of moisture induced degradation
% A : double
%   The pre-exponential constant of discoloration
% C : double
%   The constant indicating how fast the module degradates.
% Time: double
%   The time for which the degradation needs to be simulated
%
% Returns
% -------
% k_pred: double
%   The predicted degradation rate due to discoloration
%
% Developed by by Youri Blom
k = TOOLBOX_input.constants.k;
q = TOOLBOX_input.constants.q;

%Determine RMC
if TOOLBOX_input.Degradation.loadCOMSOL
    RMCFilePath = [TOOLBOX_input.Degradation.COMSOL_file,'.xlsx'];
    RMC_data = readcell(RMCFilePath,"Sheet","MoistureIngress","Range","A2:B7302");
    Time_sim = cell2mat(RMC_data(:,1));
    RMC_sim = cell2mat(RMC_data(:,2));

    %The RMC needs to be repeated for the operating after 20 years
    if Time(end) > Time_sim(end)
        N_repeat = (Time(end)-Time_sim(end))/8760;
        t_1 = find(Time_sim==Time_sim(end)-8760)+1;
        t_2 = find(Time_sim==Time_sim(end));
        delta_t_sim = Time_sim(2)-Time_sim(1);
        RMC_repeat = RMC_sim(t_1:t_2);

        Time_sim = [Time_sim; (Time_sim(end)+delta_t_sim:delta_t_sim:Time(end))'];
        RMC_sim = [RMC_sim;repmat(RMC_repeat,N_repeat,1)];

    end

    RMC = max(interp1(Time_sim,RMC_sim,Time),0);
else


end

%Peck model
k_pred = A*exp(-q*Ea/k./T(Time)).*RMC.^C;
end