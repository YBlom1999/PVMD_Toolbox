function [reflectivity,eff_Albedo] = getMaterialSpectralReflectivity(material,newLambda,libpath)
reflectivity = nan(length(material),length(newLambda),'single');
eff_Albedo = nan(length(material),'single');
%Negative values are used for specular materials
MIRROR_MAT = -1;
GLASS_MAT = -2;

amdata=csvread('am1_5_spectrum.csv');
lambdaAM15 = amdata(:,1);
specIrr15 = amdata(:,2);

for k=1:length(material)
    switch(material{k})
        case {'solar_cell','cell'}
            ref = zeros(1,length(newLambda),'single');
        case {'camera'}
            ref = zeros(1,length(newLambda),'single');
        case 'mirror'
            ref = MIRROR_MAT*ones(1,length(newLambda),'single');
        case {'glass'}
            ref = GLASS_MAT*ones(1,length(newLambda),'single');
        case 'white_glass'
            load(fullfile(libpath,'White frosted glass'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'white_paint'
            load(fullfile(libpath,'White paint'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case {'blue_paint','blue'}
            load(fullfile(libpath,'Blue paint'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'trees'
            load(fullfile(libpath,'Walnut tree'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'grass'
            load(fullfile(libpath,'Green grass'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'concrete'
            load(fullfile(libpath,'Concrete'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'pebbles'
            load(fullfile(libpath,'Pebbles'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'black'
            ref = zeros(1,length(newLambda),'single');
        case 'aluminum'
            ref = 0.4;
        case 'ground'
            load(fullfile(libpath,'Dry soil'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'new_black_mat'
            load(fullfile(libpath,'New mat'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'old_black_mat'
            load(fullfile(libpath,'Old mat'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'cable_holder'
            load(fullfile(libpath,'Black plastic'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'aluminum_profile'
            load(fullfile(libpath,'Aluminum profile'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'white_mdf'
            load(fullfile(libpath,'White MDF'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        case 'default'
            load(fullfile(libpath,'Concrete'),'lambda','specRefl');
            ref = interp1(lambda,specRefl,newLambda);
            eff_Albedo = trapz(lambdaAM15,interp1(lambda,specRefl,lambdaAM15).*specIrr15)/trapz(lambdaAM15,specIrr15);
        otherwise
            warning(['The following material is the 3D model is not recognized: ',material{k}]);
            ref = zeros(1,length(newLambda),'single');
    end
    reflectivity(k,:) = ref;
end
end