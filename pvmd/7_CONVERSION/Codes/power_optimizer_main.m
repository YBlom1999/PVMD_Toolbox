function conversion_out = power_optimizer_main(module, inverter, ...
    cable_losses, simulation, power_optimizer)
%Main function for the power optimizer
%
% Since the power optimizer is a special type of inverter, it requires of a
% separate function because SNL model is not employed.
%
% Parameters
% ----------
% module : struct
%   Electrical parameters at module level
% inverter : struct
%   User inputs related to the inverter characteristics
% cable_losses : struct
%   User inputs related to the cable losses
% simulation : struct
%   Details of the simulation
% power_optimizer : struct
%   User inputs related to the power optimizer characteristics
%
% Returns
% -------
% conversion_out : struct
%   Main results for the conversion block related to the power optimizer


[pow_opt_data, Vmpp] = get_converter_data(power_optimizer.model, ...
    simulation.periodic, module.voltage, 'pow-opt');
[pow_opt_data_STC, Vmpp_STC] = get_converter_data(power_optimizer.model, ...
    simulation.periodic, module.voltage_STC, 'pow-opt');

% Get and clean efficiencies curves as a function of input power for each
% of the input voltages previously defined
eff_curve = cell(3,1);
for j = 1:3
    eff_curve{j} = get_efficiency_curves(pow_opt_data,j);
end

eff_curve_STC = cell(3,1);
for j = 1:3
    eff_curve_STC{j} = get_efficiency_curves(pow_opt_data_STC,j);
end

% Interpolation to retrieve efficiencies at input power
if simulation.periodic
     eff_interp = cellfun(@(x) interp1(x(:,1),x(:,2),module.power),...
         eff_curve, 'UniformOutput',false);
     eff_interp_STC = cellfun(@(x) interp1(x(:,1),x(:,2),module.power_STC),...
         eff_curve_STC, 'UniformOutput',false);
else
    eff_interp = cell(inverter.num_modules,1);
    for j = 1:inverter.num_modules
        eff_interp{j} = cell(3,1);
        eff_interp{j} = cellfun(@(x) interp1(x(:,1),x(:,2),module.power{j}),...
            eff_curve, 'UniformOutput',false);
    end

    eff_interp_STC = cell(inverter.num_modules,1);
    for j = 1:inverter.num_modules
        eff_interp_STC{j} = cell(3,1);
        eff_interp_STC{j} = cellfun(@(x) interp1(x(:,1),x(:,2),module.power_STC{j},'linear','extrap'),...
            eff_curve_STC, 'UniformOutput',false);
    end
end

% Interpolation to retrieve efficiency at input voltage
pow_opt = struct();
if simulation.periodic
    eff_pow_opt = interpolate_efficiency(module.voltage,Vmpp,eff_interp);
    pow_opt.power = module.power.*eff_pow_opt/100;
    pow_opt.current = module.current;
    pow_opt.voltage = pow_opt.power./pow_opt.current;
    pow_opt.voltage(isnan(pow_opt.voltage)) = 0;


    eff_pow_opt_STC = interpolate_efficiency(module.voltage_STC,Vmpp_STC,eff_interp_STC);
    pow_opt.power_STC = module.power_STC.*eff_pow_opt_STC/100;
    pow_opt.current_STC = module.current_STC;
    pow_opt.voltage_STC = pow_opt.power_STC./pow_opt.current_STC;
    pow_opt.voltage_STC(isnan(pow_opt.voltage_STC)) = 0;
    
    Pdc_tot = sum(pow_opt.power);
    eff_PO_tot = Pdc_tot/sum(module.power)*100;
else
    eff_pow_opt = cell(inverter.num_modules,1);
    for j = 1:inverter.num_modules
        eff_pow_opt{j} = interpolate_efficiency(module.voltage{j},Vmpp,eff_interp{j});
    end
    pow_opt.power = cellfun(@(x,y) x.*y/100,eff_pow_opt,module.power,...
        'UniformOutput',false);
    pow_opt.current = module.current;    
    pow_opt.voltage = cellfun(@(x,y) x./y, pow_opt.power, pow_opt.current,...
        'UniformOutput',false);
    pow_opt.voltage = cellfun(@(x) fillmissing(x,'constant',0), ...
        pow_opt.voltage,'UniformOutput', false);
  
    
    Pdc_tot = sum([pow_opt.power{:}]);
    eff_PO_tot = Pdc_tot/sum([module.power{:}])*100;

    eff_pow_opt_STC = cell(inverter.num_modules,1);
    for j = 1:inverter.num_modules
        eff_pow_opt_STC{j} = interpolate_efficiency(module.voltage_STC{j},Vmpp_STC,eff_interp_STC{j});
    end
    pow_opt.power_STC = cellfun(@(x,y) x.*y/100,eff_pow_opt_STC,module.power_STC,...
        'UniformOutput',false);
    pow_opt.current_STC = module.current_STC;    
    pow_opt.voltage_STC = cellfun(@(x,y) x./y, pow_opt.power_STC, pow_opt.current_STC,...
        'UniformOutput',false);
    pow_opt.voltage_STC = cellfun(@(x) fillmissing(x,'constant',0), ...
        pow_opt.voltage_STC,'UniformOutput', false);
end
fprintf('The predicted power optimizer output energy yield is %0.2f Wh\n',...
    sum(Pdc_tot));
fprintf('The overall power optimizer efficiency is %0.2f %%\n',eff_PO_tot);

% Add central inverter
if power_optimizer.central_inverter
    [system,string] = compute_pv_power(pow_opt,inverter,simulation.periodic);
    [inverter_data, ~] = get_converter_data(inverter.model,nan,nan,'inv');
    
    % Retrieve efficiency data
    inv_ac_power = inverter_data(:,1);
    inv_eff = inverter_data(:,2);
    inv_dc_power = inv_ac_power./inv_eff*100;
    inv_curve = [inv_eff inv_dc_power];
    inv_curve = interpolate_first_point(inv_curve);
    
    % Fixed voltage calculations
    if power_optimizer.fixed_voltage && simulation.periodic
        Vnom_inv = inverter_data(1,3);
        pow_opt.current = pow_opt.power/Vnom_inv;
        
        string.current = pow_opt.current;
        string.voltage = string.power/string.current;
        string.voltage(isnan(string.voltage)) = 0;
        
        system.current = pow_opt.current*inverter.parallel;
        system.voltage = system.power/system.current;
        system.voltage(isnan(system.voltage)) = 0;
    end
    
    % Cable losses
    system.power = include_cable_losses(system.power,string.current,...
        system.current,cable_losses,inverter.type);

    % Inverter clipping
    system.power(system.power>inv_curve(end,2)) = inv_curve(end,2);
    
    % Interpolate inverter efficiency
    inv_eff = interp1(inv_curve(:,2),inv_curve(:,1),system.power);
    Pac = inv_eff/100.*system.power;
    Pac_tot = sum(Pac);
    inv_eff_tot = Pac_tot/sum(system.power)*100;
    
    fprintf('The predicted inverter output energy yield is %0.2f Wh [AC]\n',...
        Pac_tot);
    fprintf('The overall inverter efficiency is %0.2f %%\n',inv_eff_tot);
else
    Pac = nan; system.power = nan; inverter.model = nan;
end

if inverter.plot_graphs
    plot_pow_opt_graphs(...
        module.power, pow_opt.power,...
        eff_pow_opt,...
        eff_curve, Vmpp,...
        Pac, system.power,...
        power_optimizer.model, inverter.model, simulation.periodic)
end

conversion_out = struct();
conversion_out.Pdc = module.power;
conversion_out.Pdc_pow_opt = pow_opt.power;
conversion_out.eff_pow_opt = eff_pow_opt;
if power_optimizer.central_inverter
    conversion_out.Pac_total = Pac_tot;
    conversion_out.Pac = Pac;
    conversion_out.eff = inv_eff;
    conversion_out.eff_overall = inv_eff_tot;
end

end


function [converter_data, Vmpp] = get_converter_data(conv_model, ...
    periodic, Vdc, type)
% Get necessary data from the converter (either inverter or power
% optimizer) once the user choice has been confirmed to be suitable
%
% Parameters
% ----------
% conv_model : char
%   User-selected model for the converter
% periodic : logical
%   Whether the simulations are periodic or not
% Vdc : double/cell
%   DC voltage at module level
% type : char
%   Whether the converter is an inverter 'inv' or power optimizer 'pow-opt'
%
% Returns
% -------
% converter_data : double
%   Important data from the converter
% Vmpp : double
%   More relevant data from the converter

data_folder = get_data_path(type);

file_name = strcat(conv_model, '.xlsx');
if strcmp(type,'pow-opt') && startsWith(conv_model,'P-')
    file = fullfile(data_folder,'Solar Edge',file_name);
elseif strcmp(type,'inv') && startsWith(conv_model,'SE')
    file = fullfile(data_folder,'Solar Edge',file_name);
else
    file = fullfile(data_folder,'User Specified',file_name);
end

converter_data = xlsread(file);

if strcmp(type,'pow-opt')
    % Check maximum voltage
    Vmax_in = converter_data(4,9);
    if periodic
        Vmax_system = max(Vdc);
    else
        Vmax_system = max([Vdc{:}],[],'all');
    end

    if Vmax_system > Vmax_in
        error(['Maximum module output voltage of %0.1g V exceeds vendor',...
            ' specified maximum input voltage of %0.1g V for the selected',...
            ' power optimizer model %s'], Vmax_system, Vmax_in, conv_model)
    end
    
    % Get vendor specified input voltages
    Vmpp = converter_data(1:3,9);
else
    Vmpp = nan;
end
end

function eff_curve = get_efficiency_curves(PO_data,id)
% Get and clean three efficiency curves for the power optimizer
%
% Parameters
% ----------
% PO_data : double
%   Characteristic data from the power optimizer
% id : double
%   Curve to be obtain. Can be 1, 2, 3
% 
% Returns
% -------
% eff_curve : double
%   Efficiency curve for the power optimizer

eff_curve = PO_data(:,(2*id-1):2*id);
eff_curve = rmmissing(eff_curve);
eff_curve = interpolate_first_point(eff_curve);

end

function eff_curve = interpolate_first_point(eff_curve)
% Add a new point at y = 0 if it is zero
%
% Parameters
% -------
% eff_curve : double
%   Efficiency curve
%
% Returns
% -------
% eff_curve : double
%   Efficiency curve

if eff_curve(1,1) ~= 0
    vq = interp1(eff_curve(1:4,1),eff_curve(1:4,2),...                 
        [0 ;eff_curve(1:4,1)],'linear','extrap');
    eff_curve = [[0;eff_curve(1:end,1)] [vq;eff_curve(5:end,2)]];
end

end

function eff_PO = interpolate_efficiency(Vdc_module,Vmpp,E)
% Make an interpolation of the three efficiency curves depending on the
% voltage of the system
%
% Parameters
% ----------
% Vdc_module : double
%   DC voltage at module level
% Vmpp : double
%   Characteristic data from the power optimizer
% E : struct
%   Interpolated efficiency curves of the power optimizer
%
% Returns
% -------
% eff_PO : double
%   Efficiency of the power optimizer as a function of time

eff_PO = zeros(length(Vdc_module),1);
for i = 1:length(Vdc_module) 
    if Vdc_module(i) > Vmpp(3)
        eff_PO(i) = E{3}(i);
    elseif Vdc_module(i) > Vmpp(2)
  