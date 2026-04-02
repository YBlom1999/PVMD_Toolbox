%EXAMPLE 3: simulation of TEXTURED c-Si heterojunction cell

clear Lay Int                    %clear variables

load('AFM.mat','pyramids_20um')  %load 20x20um height map
%height map of typical pyramid texture (~5um pyramid size)

%===LAYERS===
Lay(1).med = 'air';              Lay(1).thi = inf;
Lay(2).med = 'c-Si';             Lay(2).thi = 0:1:200;
Lay(3).med = 'air';              Lay(3).thi = inf;

%===INTERFACES===
Int(1).model = 'ray';            %use RAY-optics model
Int(1).Z = pyramids_20um;        Int(1).xy = [20,20];        %20x20um map

Int(1).coat(1).med = 'ITO';      Int(1).coat(1).thi = 0.080;
Int(1).coat(2).med = 'a-Si(p)';  Int(1).coat(2).thi = 0.005;
Int(1).coat(3).med = 'a-Si(i)';  Int(1).coat(3).thi = 0.005;
%---

Int(2).model = 'ray';            %use RAY-optics model
Int(2).Z = -pyramids_20um;       Int(2).xy = [20,20];        %20x20um map
%minus was added to get upside-down height map for rear

Int(2).coat(1).med = 'a-Si(i)';  Int(2).coat(1).thi = 0.005;
Int(2).coat(2).med = 'a-Si(n)';  Int(2).coat(2).thi = 0.005;
Int(2).coat(3).med = 'Ag';       Int(2).coat(3).thi = 0.300;
%===

%---override default settings (see 'gp4_settings.m')---
S.wav = 0.300:0.050:1.200;       %wavelengths [um]
S.SM = 1;                        %output SM (see example4)

[Lay,Int,out] = GENPRO4(Lay,Int,S);  %'S' is 3rd input
%simulation output is stored in 'Lay', 'Int' and 'out'