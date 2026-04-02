function [Cost,actCfin,actE,Energy]=LCOEsc_3(Csyst,GRP,tnowyear,Pdcrat,LifeT,Com,Cmod,Exp,r,Pdd,Ey,tinv)

%% DESCRIPTION
% LCOEsc_3 FUNCTION IS USED TO SIMULATE THE THIRD SCENARIO OF LCOE MODELLING 
% 3-WHEN THE MODULES ARE REPLACED BY NEW TECHNOLOGY AFTER REACHING
% LIFETIME WITH FUTURE COSTS.
%
%----INPUT----
% Csyst - TOTAL SYSTEM COST PER kWp | Pdcrat - DC RATING OF THE SYSTEM
% GRP - GROWTH RATE PER YEAR 
% PdCrat - SYSTEM RATED POWER
% LifeT - LIFETIME OF MODULE TO REACH 80% | Com - COST OF OPERATION AND
% MAINTANANCE PER kWp 
% Cmod - EACH SECTION OF SYSTEM COST - [MODULE INVERTER STRUCTURAL_BOS
% ELECTRICAL_BOS]
% Exp - EXPECTED YEAR OF MODULES TO FUNCTION
% r - DISCOUNT RATE
% Pdd - NORMALIZED POWER WITH DEGRADATION 
% Ey - INTITAL ENERGY YIELD OF THE SYSTEM
% tinv - INVERTER LIFETIME
%
%----OUTPUT----
% COST - LCOE FOR EACH YEAR CONSIDERED AS TOTAL RECOVERY YEAR DUE TO SCENARIO 3
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
% Inverter and operation&maintenance cost intialization 
    It0=Csyst*Pdcrat ;
    Ct0=0; 
    It1=zeros(ceil(tin),1); 
    Ct1=(Com*Pdcrat); 

% Calculate number of time modules need to replaced 
    t=1:ceil(tin) ;  
    d= ~logical(rem(t,tinv)); 
    It1(d==1)= Cmod(2) ;

% Future cost of PV modules used to make a mathematical model for cost calculation 
xavg=[2014;2020;2025;2030;2035;2040;2045;2050];
yavg=[995;823;724;651;583;526;479;436];
y1=[995;830;745;684;636;595;558;516];
y2=[995;824;728;660;608;552;502;458];
y3=[995;820;716;640;565;502;452;409];
if GRP==1
%     [f1,gof1]=fit(xavg,y1,'poly7','Normalize','on','Robust','Bisquare');
     [f1,~]=fit(xavg,y1,'exp1');
     c = double( ~logical(rem(t,LifeT)) );
     z = tnowyear+find(c==1);
     It1(c==1) = It1(c==1)+ (f1(z)*Pdcrat) ;
elseif GRP==2
%     [f2,gof2]=fit(xavg,y2,'poly7','Normalize','on','Robust','Bisquare');
     [f2,~]=fit(xavg,y2,'exp1');
     c = double( ~logical(rem(t,LifeT)) );
     z = tnowyear+find(c==1);
     It1(c==1) = It1(c==1)+ (f2(z)*Pdcrat) ;
elseif GRP==3
%     [f3,gof3]=fit(xavg,y3,'poly7','Normalize','on','Robust','Bisquare');
     [f3,~]=fit(xavg,y3,'exp1');
     c = double( ~logical(rem(t,LifeT)) );
     z = tnowyear+find(c==1);
     It1(c==1) = It1(c==1)+ (f3(z)*Pdcrat) ; 
elseif GRP==4
%     [f4,gof4]=fit(xavg,yavg,'poly7','Normalize','on','Robust','Bisquare');
     [f4,~]=fit(xavg,yavg,'exp1') ;
     c = double( ~logical(rem(t,LifeT)) );
     z = tnowyear+find(c==1);
     It1(c==1) = It1(c==1)+ (f4(z)*Pdcrat) ;
end


% Discounted and Net present values of cost and energy    
    actC= (It1+Ct1)./(((1+(r/100)).^t)') ; 
    Et= (Ey*Pdcrat.*Pd); 
    Energy = sum(Et); 
    actE= Et./(((1+(r/100)).^t)') ; 
    actCfin=sum(actC)+It0; 
    actEfin=sum(actE); 
    Cost(tin,1)= actCfin./actEfin ;
    
%% 
 
%     trem=tin+1:length(Pd1); 
%     Etrem = (Ey*Pdcrat.*Pdrem) ;
%     Energyrem = sum(Etrem); 
%     actErem= Etrem./(((1+(r/100)).^trem)') ; 
%     actCfin2=sum(actC)+It0; 
%     actEfin2=sum(actE)-sum(actErem); 
%     Costsave(tin,1)= actCfin2./actEfin2 ;
%     Costsave(isempty(Costsave))=0 ;  Energyr