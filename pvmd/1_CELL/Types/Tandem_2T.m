function [Lay,Int,rback,TYPE,absmat,SubMod_ind] = Tandem_2T
%2T cell from HZB 2T with 29.15% efficiency from
%A. Al-Ashouri et al., “Monolithic perovskite/silicon tandem solar cell 
%with >29% efficiency by enhanced hole extraction,” Science, vol. 370, 
%no. 6522, pp. 1300–1309, Dec. 2020, doi: 10.1126/science.abd4016.
%Fitted by and adapted into module by MRV
      
clear Lay Int                    %clear variables
      
load('AFM.mat','pyramids_20um')  %load texture hight map from file [um]
%height map of typical pyramid texture (~5um pyramid size)
      
%===LAYERS===
Lay(1).med = 'air';                                 Lay(1).thi = inf;
Lay(2).med = 'Glass-Fe10ppmM1';                     Lay(2).thi = 3200;
Lay(3).med = 'Polyolefin-UVT';                      Lay(3).thi = 450;
%For 4T a buffer layer could be added here
Lay(4).med = 'c-Si-2015';                           Lay(4).thi = 160;
Lay(5).med = 'air';                                 Lay(5).thi = inf;
      
%===INTERFACES===
Int(1).coat(1).med = 'AF2400';               Int(1).coat(1).thi = 0.093;
Int(1).coat(2).med = 'GlassARC';             Int(1).coat(2).thi = 0.053;
      
      
Int(3).coat(1).med = 'IZO';                  Int(3).coat(1).thi = 0.085;
Int(3).coat(2).med = 'SnO2-Man18';           Int(3).coat(2).thi = 0.005;
Int(3).coat(3).med = 'C60';                  Int(3).coat(3).thi = 0.007;
Int(3).coat(4).med = 'perovskite-1.68eV-Man18-Eg-705-715';         Int(3).coat(4).thi = 0.550;
Int(3).coat(5).med = 'PTAA';                 Int(3).coat(5).thi = 0.023;
%For 4T a buffer layer could be added here
Int(3).coat(6).med = 'IWO_g';                Int(3).coat(6).thi = 0.063;
Int(3).coat(7).med = 'n-SiOx-ncSi';          Int(3).coat(7).thi = 0.111;
Int(3).coat(8).med = 'a-Si(i)';              Int(3).coat(8).thi = 0.009;
      
  %---
Int(4).coat(1).med = 'a-Si(i)';              Int(4).coat(1).thi = 0.006;
Int(4).coat(2).med = 'a-Si(p)';              Int(4).coat(2).thi = 0.012;
Int(4).coat(3).med = 'IWO_g';                  Int(4).coat(3).thi = 0.055;
Int(4).coat(4).med = 'Ag';                   Int(4).coat(4).thi = 0.3;
      
  %---
Int(3).model = 'ray';                        %use RAY-optics model
Int(3).Z = -pyramids_20um;                   Int(3).xy = [20,20];        %20x20um map
    
Int(4).model = 'ray';                        %use RAY-optics model
Int(4).Z = -pyramids_20um;                   Int(4).xy = [20,20];        %20x20um map
%minus was added to get upside-down height map for rear
  %===
      
  %===
%The reflectance of the module back is assumed to  be 0.7
rback = 0.90;
TYPE='Tan';
      
%Absorbermaterials are defined here
%In case of a tandem cell the top absorber needs to be named first
absmat = {'perovskite-1.68eV-Man18-Eg-705-715','c-Si-2015'};  
SubMod_ind = [1,1];
 end
