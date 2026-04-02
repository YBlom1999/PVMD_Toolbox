function [out] = example5()
%Current match perovskite/Si tandem by varying perovskite thickness.

%This example is a 'function' to nest a current matching function inside.
%It does not require input, so can be run in the same way as a script (F5)
%Use output '[ ] =' to pass data to Matlab workspace

load('AFM.mat','pyramids_20um')  %load 20x20um height map

%===LAYERS===
Lay(1).med = 'air';                     Lay(1).thi = inf;
Lay(2).med = 'c-Si';                    Lay(2).thi = 200;
Lay(3).med = 'air';                     Lay(3).thi = inf;

%===INTERFACES===
%interface 1: between layer 1 and 2 (air/c-Si)  
Int(1).coat(1).med = 'MgF2';            Int(1).coat(1).thi = 0.100; 
Int(1).coat(2).med = 'ITO';             Int(1).coat(2).thi = 0.040;  
Int(1).coat(3).med = 'TiO2';            Int(1).coat(3).thi = 0.040;  
Int(1).coat(4).med = 'perovskite';      Int(1).coat(4).thi = 0.400;%<-GUESS  
Int(1).coat(5).med = 'NiO';             Int(1).coat(5).thi = 0.040;  
Int(1).coat(6).med = 'a-SiOx';          Int(1).coat(6).thi = 0.040;  
Int(1).coat(7).med = 'a-Si(n)';         Int(1).coat(7).thi = 0.005;  
Int(1).coat(8).med = 'a-Si(i)';         Int(1).coat(8).thi = 0.005;  

%interface 2: between layer 2 and 3 (c-Si/air)
Int(2).model = 'ray';            %use RAY-optics model
Int(2).Z = -pyramids_20um;       Int(2).xy = [20,20];        %20x20um map          

Int(2).coat(1).med = 'a-Si(i)';         Int(2).coat(1).thi = 0.005;
Int(2).coat(2).med = 'a-Si(p)';         Int(2).coat(2).thi = 0.005;
Int(2).coat(3).med = 'Ag';              Int(2).coat(3).thi = 0.300;

%===

[tm] = fzero(@match,Int(1).coat(4).thi);   %finds zero-mismatch thickness

Int(1).coat(4).thi = tm;                   %set thickness

[Lay,Int,out] = GENPRO4(Lay,Int);          %run once more for final result


%==========================================================================
    
    function mm = match(t)
        %Set perovskite thickness, run GP4 and get mismatch.
        
        %Current matching function 'nested' in parent function (not script)
        
        Int(1).coat(4).thi = t;               %set perovskite thickness
        
        %---RUN GP4 (with SM 'trick' for speed)---
        S.SM = 1;         %output scatter matrices (both interface 1 and 2)
        [Lay,Int,~] = GENPRO4(Lay,Int,S);     %run simulation
        Int(1).SM = [];                       %delete SM interface 1
        %Note: perovskite thickness changes so SM_1 must be recalculated,
        %rear interface is unchanged so keep SM_2 and it will be reused.
        
        J1 = Int(1).coat(4).cur;              %top cell current [mA/cm2]
        J2 = Lay(2).cur;                      %bot cell current [mA/cm2]

        disp(['t=',num2str(1000*t,'%5.0f'),' nm, ',...
            'J1=',num2str(J1,'%5.2f'),' mA/cm^2, ',...
            'J2=',num2str(J2,'%5.2f'),' mA/cm^2']) %display info
       
        mm = J1-J2;                           %current mismatch
                                              %+/- bot/top current limiting
    end

end