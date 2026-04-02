y=csvread('am1_5_spectrum.csv');
lambdaAM15 = y(:,1);
specIrr15 = y(:,2);

fileList = dir('*\Irradiance Models\Spectral reflectivity library\*.mat');
fig = figure('Units','centimeters');
fig.Position(3:4) = [9,5];
ax = axes(fig);
ax.Box = 'on';
ax.NextPlot = 'add';
colors = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
lineStyle = {'-','--',':'};
nStyles = length(lineStyle);
for k=1:length(fileList)
    load(fileList(k).name,'lambda','specRefl');
    specReflAux = interp1(lambda,specRefl,lambdaAM15,'spline',0);
    albedo = trapz(lambdaAM15,specReflAux.*specIrr15)/trapz(lambdaAM15,specIrr15);
    iStyle = rem(k,nStyles);
    if iStyle == 0
        iStyle = 3;
    end
    curveName = sprintf('%s  (\\rho=%.2f)',fileList(k).name(1:end-4),albedo);
    plot(ax,lambda,100*specRefl,lineStyle{iStyle},'Color',colors(1+floor(k/nStyles-0.1),:),'LineWidth',1.1,'DisplayName',curveName);
end

legend(ax);
ax.XLim =[300 2500];
ax.XTick = 300:200:2500;
ax.XLabel.String = 'Wavelength (nm)';
ax.YLabel.String = 'Reflectivity (%)';
ax.XLabel.FontWeight = 'Bold';
ax.YLabel.FontWeight = 'Bold';
ax.FontName = 'Times New Roman';
ax.FontSize = 8;

%%