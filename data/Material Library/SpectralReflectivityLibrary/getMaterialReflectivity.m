%% Material reflectivity list
function [reflectivity,material_found] = getMaterialReflectivity(material,libpath)
reflectivity = nan(length(material),1,'single');
material_found = nan(length(material),1,'single');
%Negative values are used for specular materials
MIRROR_MAT = -1;
GLASS_MAT = -2;

amdata=csvread('am1_5_spectrum.csv');
lambdaAM15 = amdata(:,1);
specIrr15 = amdata(:,2);

for k=1:length(material)
    material_found(k) = 1;
    switch(material{k})
        case {'solar_cell','cell'}
            aux = 0;
        case {'camera'}
            aux = 0;
        case 'mirror'
            aux = MIRROR_MAT;
        case {'glass'}
            aux = GLASS_MAT;
        case 'white_glass'
            load(fullfile(libpath,'White frosted glass'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'white_paint'
            load(fullfile(libpath,'White paint'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case {'blue_paint','blue'}
            load(fullfile(libpath,'Blue paint'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'trees'
            load(fullfile(libpath,'Walnut tree'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'grass'
            load(fullfile(libpath,'Green grass'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'concrete'
            load(fullfile(libpath,'Concrete'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'pebbles'
            load(fullfile(libpath,'Pebbles'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'black'
            aux = 0;
        case 'aluminum'
            load(fullfile(libpath,'Aluminum profile'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'ground'
            load(fullfile(libpath,'Dry soil'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'new_black_mat'
            load(fullfile(libpath,'New mat'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'old_black_mat'
            load(fullfile(libpath,'Old mat'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'cable_holder'
            load(fullfile(libpath,'Black plastic'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'default'
            load(fullfile(libpath,'Concrete'),'lambda','specRefl');
            aux = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        otherwise
            aux = 0;
            material_found(k) = 0;
    end
    reflectivity(k) = aux;
end