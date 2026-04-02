function [Cost,actCfin,actE,Energy]=LCOEsc_2(Csyst,Pdcrat,LifeT,Com,Cmod,Exp,r,Pdd,Ey,tinv)
%% DESCRIPTION
% LCOESC_2 FUNCTION IS USED TO SIMULATE THE SECOND SCENARIO OF LCOE MODELLING 
% 2-WHEN THE MODULES ARE REPLACED AT SAME COST AT INSTALLATION (WITHOUT ANY
% CHANGE IN TECHNOLOGICAL) AFTER THE LIFETIME IS REACHED
%
%----INPUT----
% Csyst - TOTAL SYSTEM COST PER kWp | Pdcrat - DC RATING OF THE SYSTEM
% PdCrat - SYSTEM RATED POWER
% LifeT - LIFETIME OF MODULE TO REACH 80% | Com - COST OF OPERATION AND
% MAINTANANCE PER kWp 
% Cmod - EACH SECTION OF SYSTEM COST - [MODULE INVERTER STRUCTURAL_BOS
% ELECTRICAL_BOS]
% Exp - EXPECTED YEAR OF MODULES TO FUNCTION
% r - DISCOUNT RATE
% PdD - NORMALIZED POWER WITH DEGRADATION 
% Ey - INTITAL ENERGY YIELD OF THE SYSTEM
% tinv - INVERTER LIFETIME
%
%----OUTPUT----
% COST - LCOE FOR EACH YEAR CONSIDERED AS TOTAL RECOVERY YEAR DUE TO SCENARIO 2
% actCfin - TOTAL SYSTEM COST TILL EXPECTED PERIOD  
% actE - DISCOUNTED ENERGY TO VERIFY SYSTEM COST WITH LCOE 
% Energy - TOTAL ENERGY PRODUCED BY SYSTEM 
%% 

for tin=1:ceil(Exp)
    
% rearrange normalized power loss per year into required array size 
    Pd1=repmat(Pdd,ceil(tin/LifeT),1) ;
    Pdrem=Pd1(tin+1:end,1); 
    Pd=Pd1(1:tin,1);
    
% system cost    
    It0=Csyst*Pdcrat ;
    Ct0=0; 

% Inverter and operation&maintenance cost intialization 
    It1=zeros(ceil(tin),1); 
    Ct1=(Com*Pdcrat); 

% Calculate number of time inverter needs to be replaced 
    t=1:ceil(tin) ;  
    d= ~logical(rem(t,tinv)); 
    It1(d==1)= Cmod(2) ;
% Calculate number of time modules need to replaced 
    c = ~logical(rem(t,LifeT)); 
    It1(c==1) = It1(c==1)+ Csyst*Pdcrat ; 
    
% Discounted and Net present values of cost and energy 
    actC= (It1+Ct1)./(((1+(r/100)).^t)') ; 
    Et= (Ey*Pdcrat.*Pd); 
    Energy = sum(Et); 
    actE= Et./(((1+(r/100)).^t)') ; 
    actCfin=sum(actC)+It0; 
    actEfin=sum(actE); 
    Cost(tin,1)= actCfin./actEfin ;
    
%% 

%     Etrem = (Ey*Pdcrat.*Pdrem) ;
%     Energyrem = sum(Etrem); 
%     trem=tin+1:length(Pd1); 
% 
%     actErem= Etrem./(((1+(r/100)).^trem)') ; 
% 
%     actCfin2=sum(actC)+It0; 
%     actEfin2=sum(actE)-sum(actErem); 
% 
%     Costsave(tin,1)= actCfin2./actEfin2 ;
% 
%     Costsave(isempty(Costsave))=0 ;  Energyrem(isempty(Energyrem))=0;
%   