function [system, string] = compute_pv_power(module,inverter,periodic)
%Compute the system PV power considering the distribution of modules
%
% Parameters
% ----------
% module : struct
%   Electrical parameters at module level
% inverter : struct
%   User inputs related to the inverter characteristics
% simulation : struct
%   Details of the simulation
%
% Returns
% -------
% system : struct
%   power, voltage and current at system level (central, power optimizer
%   inverters)
% string : struct
%   power, voltage and current at string level (string, central, power
%   optimizer inverters)

system = struct();
string = struct();
if periodic
    if contains(['CEN','POW-OPT','STR'],inverter.type)
        string.current = module.current;
        string.voltage = module.voltage*inverter.series;
        string.power = module.power*inverter.series;
        string.current_STC = module.current_STC;
        string.voltage_STC = module.voltage_STC*inverter.series;
        string.power_STC = module.power_STC*inverter.series;
    end
    if contains(['CEN','POW-OPT'],inverter.type)
        system.power = module.power*inverter.parallel*inverter.series;
        system.current = module.current*inverter.parallel;
        system.voltage = module.voltage*inverter.series;
        system.power_STC = module.power_STC*inverter.parallel*inverter.series;
        system.current_STC = module.current_STC*inverter.parallel;
        system.voltage_STC = module.voltage_STC*inverter.series;
    end
else
    module = structfun(@cell2mat,module,'UniformOutput',false);
    resize_shape = [size(module.power,1)/inverter.num_modules,...
        inverter.num_modules];
    module_operating = {};
    module_operating.power = module.power;
    module_operating.voltage = module.voltage;
    module_operating.current = module.current;
    system = structfun(@(x) reshape(x,resize_shape),module_operating,...
        'UniformOutput',false);
    system.power_STC = module.power_STC';
    system.current_STC = module.current_STC';
    system.voltage_STC = module.voltage_STC';
    if contains(['STR','CEN','POW-OPT'],inverter.type)
        comb = combinations(inverter.num_strings,inverter.num_modules);
        string.power = cell(inverter.num_strings,1);
        string.voltage = cell(inverter.num_strings,1);
        string.current = cell(inverter.num_strings,1);
        for i = 1:inverter.num_strings
            comb_str = comb(i,:);
            [string.power{i}, string.voltage{i}] = mismatch(...
                system.voltage,...
                system.current,...
                comb_str);
        end
        string.current = cellfun(@(x,y) x./y,...
            string.power,...
            string.voltage,...
            'UniformOutput',false);
        string.current = cellfun(@(x) fillmissing(x,'constant',0), ...
            string.current,'UniformOutput', false);
        
        string.power_STC = cell(inverter.num_strings,1);
        string.voltage_STC = cell(inverter.num_strings,1);
        string.current_STC = cell(inverter.num_strings,1);
        for i = 1:inverter.num_strings
            comb_str = comb(i,:);
            [string.power_STC{i}, string.voltage_STC{i}] = mismatch(...
                system.voltage_STC,...
                system.current_STC,...
                comb_str);
        end
        string.current_STC = cellfun(@(x,y) x./y,...
            string.power_STC,...
            string.voltage_STC,...
            'UniformOutput',false);
        string.current_STC = cellfun(@(x) fillmissing(x,'constant',0), ...
            string.current_STC,'UniformOutput', false);


        if strcmp(inverter.type,'STR')
            system = struct();
        else
            % Calculate voltage mismatch from parallel strings
            system = structfun(@(x) sum([x{:}],2), string,'UniformOutput',false);
            system.voltage = system.power./system.current;
            system.voltage(isnan(system.voltage)) = 0;
        end
    end
end

end


function comb = combinations(sz,nummod)
% Output indices of panels present in a string
%
% Parameters
% ----------
% sz : double
%   number of strings in the system
% nummod :double
%   number of modules in the system
%
% Returns
% -------
% comb : double
%   panel indices of a single string for string selection. Each row has the
%   indices of panels in a string
%
% Author: K Ganapathi Subramanian

if round(nummod/sz) ~= nummod/sz
   error('Division to %d strings not possible!',sz); 
elseif sz == 1
   comb = 1:nummod;
elseif sz == nummod
    %micro-inverter case
    comb = (1:nummod)';
else
    s = 1:nummod;
    comb = (reshape(s,nummod/sz,sz))';
end
end


function [Pdc_string,Vdc_string] = mismatch(Vdc_mpp,Idc_mpp,comb)
%Calculate the mismatch losses from the PV system
% Squared approx for MPP of a panel. The string current is the same as that
% of the module
% 
% Parameters
% ----------
% Vdc_mpp : double
%   MPP voltages of each panel in the string
% Idc_mpp : double
%   MPP currents of each panel in the string
% comb : double
%   panel indices of a single string for string selection. Each row has the
%   indices of panels in a string
% 
% Returns
% -------
% Pdc_string : double
%   string DC power
% Vdc_string : double
%   string DC voltage
%
% Author: K Ganapathi Subramanian

Vdc_mpp = Vdc_mpp(:,comb);
Idc_mpp = Idc_mpp(:,comb);

% Sort string panel currents
[I,idx1] = sort(Idc_mpp,2,'descend');
V = zeros(size(Vdc_mpp));
% Sort voltages as per sorted current index
for k = 1:size(Vdc_mpp,1)
    for l = 1:size(Vdc_mpp,2)
        m = idx1(k,l);
        V(k,l) = Vdc_mpp(k,m);
    end
end
% Find cumulative string voltage and power from each module
V = cumsum(V,2); 
P = I.*V;

% Find module producing maximum power
[~,idx2] = max(P,[],2); 
Pdc_string = zeros(size(P,1),1);
Vdc_string = zeros(size(P,1),1);
for time_i = 1:size(P,1)
    Idc_string = I(time_i,idx2(time_i));
    Pdc_string(time_i,1) = P(time_i,idx2(time_i));
    Vdc_string(time_i,1) = Pdc_string(time_i,1)/Idc_string;
end
Vdc_string(isnan(Vdc_string)) = 0;
end