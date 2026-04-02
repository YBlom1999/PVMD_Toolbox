function aplot(out)
% Make area plot based on GP4 simulation output. Stacking order can be
% modified to 'absorber layer(s) at the bottom' (for comparison with EQE).
% User can edit material list in 'colorful'(see bottom of this file)

%---INPUT---
%out            structure with GP4 simulation output (in logical order)
%out.leg        layer/coating names(for legend)
%out.abp        absorptance in every (slice of) layer/coating [-]
%out.cur        corresponding implied photocurrent [mA/cm2]
%out.wav        wavelength range [um]

L = length(out.ix);                 %nr of layers/coatings
W = length(out.wav);                %nr of wavelengths
data = zeros(L,W+5);                %data matrix (abs vs wl, current, rgb, a) 

out.leg{1} = 'reflected';           %change first medium 'air' to 'reflected'
out.leg{L} = 'transmitted';         %change final medium to 'transmitted'
for l = 1:L                         %for every layer
    i = L-l+1;                      %plotting index (reverse order for area plot!)
    data(i,1:W) = sum(out.abp(out.ix{l},:),1);  %sum absorptance vs wavelength over all slices
    cur = sum(out.cur(out.ix{l}));  %sum current over all slices [mA/cm2]
    data(i,W+1) = cur;              %add to matrix
    [c,a] = colorful(out.leg{l});   %lookup layer plotting color and absorber status
    data(i,W+(2:5)) = [c,a];        %add to matrix color (rgb) and absorber status
    leg{i,1} = [out.leg{l},' (',num2str(cur,'%4.1f'),' mA/cm^2)'];  %#ok<AGROW> %add current to string
end           

%---move absorber layers to bottom of area plot and merge if multipe (tandem)---
f = find(data(:,W+5)==1);          %find absorber layers
if ~isempty(f)                     %if absorber layers found
    n = max(numel(f)==1,0.3);      %color modifier
    dvec = [sum(data(f,1:W+1),1),mean(data(f,W+(2:4)),1).^n,0]; %create abs sum data vector
    data = [dvec;data];            %add it to the matrix
    leg = [[leg{f}];leg];          %add legend entry for combined area (won't be visible)
end

%---make area plot ---
figure(101)                        %
clf                                %clear figure (if exists)
ixa = 1:size(data,1);              %index to all layers (incl. abs and sum)
ixa(f+1) = [];                     %remove the absorber layers
ha=area(1000*out.wav,data(ixa,1:W).');               %make area plot
set(ha,{'FaceColor'},num2cell(data(ixa,W+(2:4)),2)); %set specified area colors
il = flipud(find(data(ixa,W+1)>0.03));                 %find layers with J > 0.03 mA/cm2
hl = legend(ha(il),leg{ixa(il)});                    %legend (show only those)

%---make line plot---
if numel(f)>1                       %only if tandem
    hold on                         %plot over area plot
    hp = plot(1000*out.wav,data(f+1,1:W),'LineWidth',3); %make line plot for absorbers
    set(hp,{'color'},num2cell(data(f+1,W+(2:4)),2))      %set specified line colors
    il(end) = [];                                        %remove sum area from legend
    hl = legend([ha(il),fliplr(hp')],[leg(ixa(il));flipud(leg(f+1))]);     %legend (area + lines)
end

%---make graph look good---
set(hl,'Location','NorthEastOutside','FontName','Calibri','FontSize',12);
grid on
set(gca,'Layer','Top','FontName','Calibri','fontsize',16)
xlabel('wavelength [nm]','fontsize',18)
ylabel('absorptance [-]','fontsize',18)
xlim(1000*[out.wav(1),out.wav(end)]);
ylim([0,1]);
%--------------------------------------------------------------------------
    function [c,a] = colorful(mat_name)
        %For plotting purposes, lookup for a given material the rgb color 
        %coordinate and whether it is an absorber layer (plotted at bottom)
        
        %---list of materials (USER CAN EDIT BELOW)---
        m( 1).name = 'reflected';   m( 1).rgb = [255,255,255]/255;      m( 1).abs = 0;
        m( 2).name = 'ITO';         m( 2).rgb = [208,242,242]/255;      m( 2).abs = 0;
        m( 3).name = 'TiO2';        m( 3).rgb = [209,239, 84]/255;      m( 3).abs = 0;
        m( 4).name = 'NiO';         m( 4).rgb = [232,155, 39]/255;      m( 4).abs = 0;
        m( 5).name = 'a-Si(n)';     m( 5).rgb = [100,100,153]/255;      m( 5).abs = 0;
        m( 6).name = 'a-Si(i)';     m( 6).rgb = [126,100,126]/255;      m( 6).abs = 0;
        m( 7).name = 'a-Si(p)';     m( 7).rgb = [153,100,100]/255;      m( 7).abs = 0;
        m( 8).name = 'Ag';          m( 8).rgb = [ 51, 51, 51]/255;      m( 8).abs = 0;
        m( 9).name = 'c-Si';        m( 9).rgb = [160,160,170]/255;      m( 9).abs = 1;
        m(10).name = 'perovskite';  m(10).rgb = [191,156, 73]/255;      m(10).abs = 1;
        m(11).name = 'FTO';         m(11).rgb = [209,239, 84]/255;      m(11).abs = 0;
        m(12).name = 'PCBM';        m(12).rgb = [232,155, 39]/255;      m(12).abs = 0;
        m(13).name = 'PTAA';        m(13).rgb = [232,155, 39]/255;      m(13).abs = 0;
        m(14).name = 'a-SiOx';      m(14).rgb = [126,100,126]/255;      m(14).abs = 0;
        m(15).name = 'spiro-OMeTAD';m(15).rgb = [209,239, 84]/255;      m(15).abs = 0;
        m(16).name = 'transmitted'; m(16).rgb = [232,232,255]/255;      m(16).abs = 0;
        %Setting m().abs = 0 for all layers gives good old 'plain' stacking order.
        %---(USER CAN EDIT ABOVE)---
        
        x = find(strcmp(mat_name,{m.name}),1);  %find on list
        if isempty(x)            %if material not found
            rng(sum(mat_name*1)) %reproducable (string ascii sum as seed)
            c = rand(1,3);       %random color if material is not found
            a = 0;               %assume not absorber layer
        else                     %if material found
            c = m(x).rgb;        %take specified color
            a = m(x).abs;        %and specified absorber state
        end
    end
end