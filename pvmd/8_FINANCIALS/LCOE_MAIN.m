function [LCOE_output, TOOLBOX_input] = LCOE_MAIN(TOOLBOX_input,ELECTRIC_output,DEGRADATION_output)
% Description
%THIS FUNCTION CALCULATES LCOE IN 3 SCENARIOS
% 1-WHEN SAME MODULE IS USED FOR EXPECTED YEAR WITHOUT ANY REPLACEMENT
% 2-WHEN THE MODULES ARE REPLACED BY SAME MODULE (WITHOUT ANY CHANGE
% TECHNOLOGICAL CHANGES) AFTER THE LIFETIME IS REACHED
% 3-WHEN THE MODULES ARE REPLACED BY NEW TECHNOLOGY AFTER REACHING
% LIFETIME.

% Required inputs
% Inputs include number of modules to calculate rated power, discount rate,
% cost of system, maximum expected period, inverter lifetime and year of
%installation(required for calculation of future cost of modules)

if isfield(DEGRADATION_output,'Output')
    warning('LCOE calculations with User defined Panel Degradation do not represent the correct values!');
    Tworking = input('Define Number of Working Years [Max 50]: ');
    dg_usr = input('Define Panel Degradation Rate [%]: ');
    NORMALIZEDPOWER_output.PRtotal = (1-dg_usr/100).^(0:Tworking-1)';
    NORMALIZEDPOWER_output.Lifetime = Tworking;
    NORMALIZEDPOWER_output.PRworking_years = (1-dg_usr/100).^(0:Tworking-1)';
else
    NORMALIZEDPOWER_output = DEGRADATION_output.NORMALIZEDPOWER_output;
end
nummod = TOOLBOX_input.Conversion.Parallel_Modules*...
    TOOLBOX_input.Conversion.Parallel_Modules;
if ~TOOLBOX_input.script %GUI VERSION
    TOOLBOX_input.runFinancialPart = true;
    prompt = {'Discount ratio[%]',...
        'Module cost[Eur/kWp]','Inverter cost[Eur/kWp]',...
        'Structural BOS[Eur/kWp]','Electrical BOS[Eur/kWp]',...
        'O&M cost[Eur/kWp]:',...
        'Inverter Lifetime','Installation year'};
    dlgtitle = 'Input';
    dims = [1 50];
    definput = {'5','395.47','151.45','84.14','193.53','16.8','50','13','2021'};
    details = inputdlg(prompt,dlgtitle,dims,definput);
    Values=str2double(details);
    Values = [nummod;Values];
    
    TOOLBOX_input.FinancialPart.DiscountRate = Values(2);
    TOOLBOX_input.FinancialPart.Modulecost = Values(3);
    TOOLBOX_input.FinancialPart.Invertercost = Values(4);
    TOOLBOX_input.FinancialPart.StructuralBOScost= Values(5);
    TOOLBOX_input.FinancialPart.ElectricalBOScost = Values(6);
    TOOLBOX_input.FinancialPart.OMcost = Values(7);
    TOOLBOX_input.FinancialPart.Inverterlifetime = Values(8);
    TOOLBOX_input.FinancialPart.Installyear = Values(9);
else %SCRIPT VERSION
    Values(1) = nummod;
    Values(2) = TOOLBOX_input.FinancialPart.DiscountRate;
    Values(3) = TOOLBOX_input.FinancialPart.Modulecost;
    Values(4) = TOOLBOX_input.FinancialPart.Invertercost;
    Values(5) = TOOLBOX_input.FinancialPart.StructuralBOScost;
    Values(6) = TOOLBOX_input.FinancialPart.ElectricalBOScost;
    Values(7) = TOOLBOX_input.FinancialPart.OMcost;
    Values(8) = TOOLBOX_input.FinancialPart.Inverterlifetime;
    Values(9) = TOOLBOX_input.FinancialPart.Installyear;
end

% Select financial scenario - Scenario 1 - same module till expected
% functioning period, Scenario 2- degradated modules are replaced by new
% modules after lifetime (Scenario 2 is divided into subscenarios with and without growth rate per year)

if ~TOOLBOX_input.script
    Str={'Select financial scenario'} ;
    S={'Same module till required years','New technology modules replaced till required years'};
    result=listdlg('PromptString',Str,'ListString',S,'SelectionMode','single','ListSize',[300 70]);
    TOOLBOX_input.FinancialPart.scenario = result;
else
    result = TOOLBOX_input.FinancialPart.scenario;
end

% Assinging the values to the variables

% STC Power and energy to calculate energy density
Pdc=ELECTRIC_output.DCP; % rated power
if TOOLBOX_input.runPeriodic
    PSTC=ELECTRIC_output.P_STC*1E-03; % kW conversion
    Energy=sum(Pdc)*1E-03 ; % kWh conversion
    %Q: DO WE INCLUDE TOTAL NUMBER OF MODULES HERE?
    inter=NORMALIZEDPOWER_output.Lifetime;
    Pratio=NORMALIZEDPOWER_output.PRtotal;
    Pdworking=NORMALIZEDPOWER_output.PRworking_years;
    Pdcrat=PSTC*nummod; % system rated power
    Ey=nummod*(Energy/PSTC) ; % system energy specific yield (KWh/KWp)
else %non-periodic simulations
    PSTC = cellfun(@(x) x*1e-03,ELECTRIC_output.P_STC);
    Energy = cellfun(@sum,Pdc);
    inter= min(cellfun(@min,NORMALIZEDPOWER_output.Lifetime));
    len = min(cellfun(@length,NORMALIZEDPOWER_output.PRtotal));
    Pratio = cellfun(@(x) x(1:len),...
        NORMALIZEDPOWER_output.PRtotal,'UniformOutput',false);
    Pratio = cell2mat(Pratio);
    Pratio = reshape(Pratio,len,length(Pratio)/len);
    Pratio = mean(Pratio,2);
    Pdworking= cell2mat(NORMALIZEDPOWER_output.PRworking_years);
    Pdworking = reshape(Pdworking,length(Pdworking)/nummod,nummod);
    Pdworking = mean(Pdworking,2);
    Pdcrat= sum(PSTC); % system rated power
    Ey= sum(Energy./PSTC); % system energy specific yield (KWh/KWp)
end

%  Redefining the variables
r=Values(2);
Cmod=Values(3:6);  Csyst=sum(Values(3:6));
Com=Values(7);
Exp= size(Pdworking,1); ...
    %expected period is equal to number of years defined in DEGRADATION
tinv=Values(8);
tnowyear=Values(9);
% Financial calculations
% LCOE output represents LCOE per year as recovery period, SystemCost_dis -
% discounted over all system cost, Energy_dis - discounted energy cost to
% calculate verify system cost with LCOE, Energy - overall energy produced
% by system till expected functioning period



LifeT = ceil(inter);
Pdd = [1 ; Pratio(1:LifeT-1)]; % Intial values is 1 to evaluvate year without degradation

if result==1
    [Cost_1,actCfin_1,actE_1,Energy_1]=LCOEsc_1(Csyst,Pdcrat,LifeT,Com,Cmod,Exp,r,Pdworking,Ey,tinv);
    LCOE_output.S1.LCOE_1=Cost_1;
    LCOE_output.S1.SystemCost_dis_1=actCfin_1;
    LCOE_output.S1.Energy_dis_1=actE_1;
    LCOE_output.S1.Energy_1=Energy_1;
    if ~TOOLBOX_input.script
        LCOE_plots(result,Cost_1);
    end
    
elseif result==2
    
    % Growth rate per year selection to evaluvate future module costs - 4
    % options are presented 5%,7.5%, 10% and average (mean of 5%,7.5%,10%
    % and PV breakthrought scenario)
    if ~TOOLBOX_input.script
        dlgQuestion={'Select Growth rate per year'} ;
        rate={'5%(pessimistic)','7.5%(intermediate)','10%(optimistic)','Avg(recommended)'};
        GRP=listdlg('PromptString',dlgQuestion,'ListString',rate,'SelectionMode','single','ListSize',[300 250]);
        TOOLBOX_input.FinancialPart.GRP = GRP;
    else
        GRP = TOOLBOX_input.FinancialPart.GRP;
    end
    
    [Cost_1,actCfin_1,actE_1,Energy_1]=LCOEsc_1(Csyst,Pdcrat,LifeT,Com,Cmod,Exp,r,Pdworking,Ey,tinv);
    LCOE_output.S1.LCOE_1=Cost_1;
    LCOE_output.S1.SystemCost_dis_1=actCfin_1;
    LCOE_output.S1.Energy_dis_1=actE_1;
    LCOE_output.S1.Energy_1=Energy_1;
    
    [Cost_2,actCfin_2,actE_2,Energy_2]=LCOEsc_2(Csyst,Pdcrat,LifeT,Com,Cmod,Exp,r,Pdd,Ey,tinv);
    LCOE_output.S2.LCOE_2=Cost_2;
    LCOE_output.S2.SystemCost_dis_2=actCfin_2;
    LCOE_output.S2.Energy_dis_2=actE_2;
    LCOE_output.S2.Energy_2=Energy_2;
    
    
    [Cost_3,actCfin_3,actE_3,Energy_3]=LCOEsc_3(Csyst,GRP,tnowyear,Pdcrat,LifeT,Com,Cmod,Exp,r,Pdd,Ey,tinv);
    LCOE_output.S3.LCOE_3=Cost_3;
    LCOE_output.S3.SystemCost_dis_3=actCfin_3;
    LCOE_output.S3.Energy_dis_3=actE_3;
    LCOE_output.S3.Energy_3=Ener