function [Pmpp, Impp, Vmpp, Isc, Voc]=CreatingModuleDCOutput(V_module, I,time,Mod_type,Tracking_type,minVoltage)
% CreatingModuleDCOutput calculates the DC output power of the
% module and determines the IV characteristics.
%
% Based on the module IV curve, the maximum power point is found and the IV
% characteristics are determined
%
% Parameters
% ----------
% V_module : double
%   The voltage of the module IV curve
% I : double
%   The current of the module IV curve
% time : double
%   The time for which the DC output must be calculated
% Mod_type : double
%   The type of module (Non-Rec/ Rec/ Butterfly)
% Tracking_type : struct
%   The type of the maximum power point tracker (Global/ local)
% minVoltage : double
%   The lowest value the voltage is allowed to reach
%
% Returns
% -------
% Pmpp : double
%   The maximum power point power of each moment in time
% Impp : double
%   The maximum power point current of each moment in time
% Vmpp : double
%   The maximum power point voltage of each moment in time
% Isc : double
%   The short circuit current of each moment in time
% Voc : double
%   The open circuit voltage of each moment in time
%
% Performed by unknown (M.R Vogt, A. Nour)
% Commented by Y. Blom

if strcmp(Mod_type,'REC') %for reconfigurable modules
    [Pmpp, Impp, Vmpp, Isc, Voc] = ModuleOutputReconfigurable(V_module,I);

else %for standard (non-reconfigurable) modules
    if strcmp(Tracking_type,'Global')
        [Pmpp, Impp, Vmpp, Isc, Voc] = ModuleOutputGlobalTracker(V_module,I,time,minVoltage);
    elseif strcmp(Tracking_type,'Local')
        [Pmpp, Impp, Vmpp, Isc, Voc] = ModuleOutputLocalTracker(V_module,I,time,minVoltage);

    end


end

if length(V_module(:,1))>time %4T
    Isc=reshape(Isc,time,2);
    Impp=reshape(Impp,time,2);
    Vmpp=reshape(Vmpp,time,2);
    Voc=reshape(Voc,time,2);
end
%%Ensure all vectors are row vectors
Isc=reshape(Isc,time,numel(Isc)/time);
Impp=reshape(Impp,time,numel(Impp)/time);
Voc=reshape(Voc,time,numel(Voc)/time);
Vmpp=reshape(Vmpp,time,numel(Vmpp)/time);
Pmpp=reshape(Pmpp,time,numel(Pmpp)/time);
end

function [Pmpp, Impp, Vmpp, Isc, Voc] = ModuleOutputReconfigurable(V_module,I)
% ModuleOutputReconfigurable calculates the DC output power of the
% module and determines the IV characteristics for reconfigurable modules.
%
% Based on the module IV curve, the maximum power point is found and the IV
% characteristics are determined for reconfigurable modules
%
% Parameters
% ----------
% V_module : double
%   The voltage of the module IV curve
% I : double
%   The current of the module IV curve
%
% Returns
% -------
% Pmpp : double
%   The maximum power point power of each moment in time
% Impp : double
%   The maximum power point current of each moment in time
% Vmpp : double
%   The maximum power point voltage of each moment in time
% Isc : double
%   The short circuit current of each moment in time
% Voc : double
%   The open circuit voltage of each moment in time
%
% Written by Y. Blom
Pmpp=zeros(length(V_module(:,1)),1);
Vmpp=zeros(length(V_module(:,1)),1);
Impp=zeros(length(V_module(:,1)),1);
Isc=zeros(length(V_module(:,1)),1);
Voc=zeros(length(V_module(:,1)),1);

for i=1:length(V_module(:,1))
    P=V_module(i,:).*I;
    [v,r]=max(P);
    Pmpp(i)=v;
    Vmpp(i)=V_module(i,r);
    Impp(i)=I(r);
end

for i=1:length(V_module(:,1))
    for j=1:size(I,2)-1
        if V_module(i,j)>0 && V_module(i,j+1)==0
            Isc(i)=I(j+1);
            break
        end
    end
end

for i=1:length(V_module(:,1))
    if Isc(i)>0
        Voc(i)=interp1(I,V_module(i,:),0);
    end
end
end

function [Pmpp, Impp, Vmpp, Isc, Voc] = ModuleOutputGlobalTracker(V_module,I,time,minVoltage)
% ModuleOutputGlobalTracker calculates the DC output power of the
% module and determines the IV characteristics with global MPPT.
%
% Based on the module IV curve, the maximum power point is found and the IV
% characteristics are determined with global MPPT
%
% Parameters
% ----------
% V_module : double
%   The voltage of the module IV curve
% I : double
%   The current of the module IV curve
% time : double
%   The time for which the DC output must be calculated
% minVoltage : double
%   The lowest value the voltage is allowed to reach
%
% Returns
% -------
% Pmpp : double
%   The maximum power point power of each moment in time
% Impp : double
%   The maximum power point current of each moment in time
% Vmpp : double
%   The maximum power point voltage of each moment in time
% Isc : double
%   The short circuit current of each moment in time
% Voc : double
%   The open circuit voltage of each moment in time
%
% Written by Y. Blom
Pmpp=zeros(length(V_module(:,1)),1);
Vmpp=zeros(length(V_module(:,1)),1);
Impp=zeros(length(V_module(:,1)),1);
Isc=zeros(length(V_module(:,1)),1);
Voc=zeros(length(V_module(:,1)),1);

for i=1:length(V_module(:,1))
    %Find which voltages can be reached
    if max(V_module(i,:)) == 0; continue; end
    Reachable = find(V_module(i,:) > minVoltage);
    if isempty(Reachable); continue; end
    P=V_module(i,Reachable).*I(Reachable);
    [v,r]=max(P);
    Pmpp(i)=v;
    Vmpp(i)=V_module(i,Reachable(r));
    Impp(i)=I(Reachable(r));
end

if length(Pmpp)>time
    Pmpp=reshape(Pmpp,time,2);
    Pmpp=sum(Pmpp,2);
end

for i=1:length(V_module(:,1))
    [~,un_index] = unique(V_module(i,:));
    if length(un_index) > 1
        Isc(i) = interp1(V_module(i,un_index),I(un_index),0,'linear','extrap');
        if Isc(i)>0
            Voc(i)=interp1(I,V_module(i,:),0);
        end
    end
end

end

function [Pmpp, Impp, Vmpp, Isc, Voc] = ModuleOutputLocalTracker(V_module,I,time,minVoltage)
% ModuleOutputGlobalTracker calculates the DC output power of the
% module and determines the IV characteristics with global MPPT.
%
% Based on the module IV curve, the maximum power point is found and the IV
% characteristics are determined with global MPPT
%
% Parameters
% ----------
% V_module : double
%   The voltage of the module IV curve
% I : double
%   The current of the module IV curve
% time : double
%   The time for which the DC output must be calculated
% minVoltage : double
%   The lowest value the voltage is allowed to reach
%
% Returns
% -------
% Pmpp : double
%   The maximum power point power of each moment in time
% Impp : double
%   The maximum power point current of each moment in time
% Vmpp : double
%   The maximum power point voltage of each moment in time
% Isc : double
%   The short circuit current of each moment in time
% Voc : double
%   The open circuit voltage of each moment in time
%
% Written by Y. Blom
Pmpp=zeros(length(V_module(:,1)),1);
Vmpp=zeros(length(V_module(:,1)),1);
Impp=zeros(length(V_module(:,1)),1);
Isc=zeros(length(V_module(:,1)),1);
Voc=zeros(length(V_module(:,1)),1);

P=V_module(1,:).*I;
[v,r]=max(P);
Pmpp(1)=v;
Vmpp(1)=V_module(1,r);
Impp(1)=I(r);
% fig = figure;

for i=2:length(V_module(:,1))
    if max(V_module(i,:)) == 0; continue; end
    Reachable = find(V_module(i,:) > minVoltage);
    if isempty(Reachable); continue; end

    P=V_module(i,Reachable).*I(Reachable);
    if max(P) == 0; continue; end
    Max_found = false;
    [~,r] = min(abs(V_module(i,Reachable)-Vmpp(i-1)));
    while ~Max_found
        %If the current point is the higher than its neigbhour points, the
        %current point is selected as maximum
        if P(r) >= P(max(r-1,1)) && P(r) >= P(min(r+1,length(P)))
            Pmpp(i) = P(r);
            Vmpp(i) = V_module(i,Reachable(r));
            Impp(i) = I(Reachable(r));
            Max_found =  true;
        %If the point to the left is higher, the MPPT moves towards the
        %left
        elseif P(r) < P(max(r-1,1))
            r = max(r-1,1);
        %If the point to the right is higher, the MPPT moves towards the
        %right
        elseif P(r) < P(min(r+1,length(P)))
            r = min(r+1,length(P));
        end
    end

%     clf(fig);
%     hold on; box on;
%     plot(V_module(i,:),P)
%     plot(Vmpp(i),Pmpp(i),'x','MarkerSize',15)
%     xlabel('Voltage [V]')
%     ylabel('Power [W]')
%     xlim([0, 55])
%     ylim([0, 1.2*max(P)]);
%     title(append('t = ',num2str(i/6)));
%     drawnow
%     pause(0.05)
end

if length(Pmpp)>time
    Pmpp=reshape(Pmpp,time,2);
    Pmpp=sum(Pmpp,2);
end

for i=1:length(V_module(:,1))
    Isc_ind = find(V_module(i,:) == 0,1);
    Isc(i) = I(Isc_ind);
    if Isc(i)>0
        Voc(i)=interp1(I,V_module(i,:),0);
    end
end

end