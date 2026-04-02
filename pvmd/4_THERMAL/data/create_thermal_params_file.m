function create_thermal_params_file()
% Height of the anenometer [m]
anenometer_height = 10;
% Thermal conductivity of glass [W/(m K)]
glass_conduct = 1;
% Emissivity of glass [-]
glass_emissivity = 0.89;
% Stefan-Boltzmann constant [W/(m2 K4)]
boltzmann_ct =  5.6697*1e-8;

save('thermal_params.mat','anenometer_height','glass_conduct',...
    'glass_emissivity','boltzmann_ct')
end