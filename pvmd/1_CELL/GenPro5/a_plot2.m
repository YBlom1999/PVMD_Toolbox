function a_plot2(tab,ax,varargin)
% Make area plot based on GENRPO simulation output. Absorber layer is
% plotted at bottom, and merged if there are multiple (e.g. tandem).

% This is the default plot in GENPRO. Will also appear in app (Scenario A).
% Plots for other scenarios are created by GENPRO_shell.

%---INPUT---
% tab       table with a row for each slice and column for each wavelength
% ax        axis where figure will be plotted (0 = don't plot)

if isnumeric(ax)        %if ax is not axis handle (e.g. provide by GUI)
    if ax==0
        return          %silent mode --> don't plot
    else                %any positive number
        figure(ax)
        ax = gca;       %use other axis (in existing or new figure) 
    end
end

cla(ax,"reset");            %clear axis (removes any previous plots)

%---format table---
tab(startsWith(tab{:,1},'*'),:) = [];   % remove texture rows
tab = unslice(tab,[2,4:size(tab,2)]);   % merge individual slices

% get wavelength range from table header
wavelength_nm = str2double(extractBetween(tab.Properties.VariableNames(5:end),'=','nm'));
%wavelength_nm = str2double(tab.Properties.VariableNames(6:end));

R = size(tab,1);                        % nr of rows in table

[rgb,abs_file] = read_rgb(tab.Layer); % read rgb & a from file\
% note: done here (so open xlsx file only if plotting)
% note: colors are read from settings file twice. Once for creating input 
% plot (device cross section), and again here for creating output plot. 
% Reading file twice takes a bit of time, but this is the most robust way.

tab = [tab,table(rgb)];                 % add rgb column (after final wl)

%add column with absorber status (after rgb column)
tab = [tab,table(abs_file)];            % lookup absorbers in settings file
 
tab = sortrows(tab,size(tab,2));        % sort by final column 
                                        % (so absorbers move down)

ia = find(tab{:,end}==1);               % find absorbers

if numel(ia)>1                          % if more than one absorber

    %absorber areas are 'merged'. Lines indicate individual contribution.
    tab(R+1,:) = tab(R,:);              % add extra row
    tab(R+1,"Layer") = {'absorbers:'};  % assign material name

    tab{R+1,5:end-2} = sum(table2array(tab(ia,5:end-2)),1);  % sum abs.
    tab{R+1,"J [mA/cm2]"}     = sum(table2array(tab(ia,"J [mA/cm2]")),1); % sum J
    tab{R+1,"rgb"}   = color_mixer(table2cell(tab(ia,"rgb")));   % mix rgb

    lin = tab(ia,:);                    % take indivial absorber lines
    leg2 = strcat(lin.Layer,{' ('},num2str(lin.("J [mA/cm2]"),'%4.1f'),{' mA/cm2)'});

    tab(ia,:) = [];                     % remove lines from original table
end

R = size(tab,1);                        % update nr of rows in table
AA = table2array(tab(:,5:end-2));       % get absorptance data matrix
%note: absorptance data is in columns 6 till 1 before last.

%leg  = strcat(tab.Layer,{' ('},num2str(tab(:,5),'%4.1f'),{' mA/cm2)'});
leg  = strcat(tab.Layer,{' ('},num2str(tab{:,4},'%4.1f'),{' mA/cm2)'});
%TODO: use tab.("J[mA/cm2]") notation, which is more robust

il = find(tab{:,4}>0.03);                  % is large (current > 0.03 mA/cm2)
il2 = R+1-il;                           % reversed (counting from bottom)

%---make area plot ---
ha = area(ax,wavelength_nm,flipud(AA).'); %make area plot
set(ha,{'FaceColor'},flipud(tab.rgb));    %set area colors

%---add legend---
hl = legend(ax,ha(il2),leg(il));             %legend of area (not lines)
%ALL layers are ploted, legend only shows layers for which current is large

if numel(ia)>1                            %if multiple absorbers
    %---add line plot---
    hold(ax,"on")
    hp = plot(ax,wavelength_nm,lin{:,5:end-2});    %plot absorber lines
    set(hp,{'color'},lin.rgb,'LineWidth',3);    %set color & width
    hl = legend(ax,[ha(il2),hp'],[leg(il);leg2]); %update legend
end

%---make graph look good---
set(hl,'Location','NorthEastOutside','FontName','Calibri','FontSize',12);
grid(ax,'on')
set(ax,'Layer','Top','FontName','Calibri','fontsize',16)
xlabel(ax,'wavelength [nm]','fontsize',18)
ylabel(ax,'absorptance [-]','fontsize',18)
xlim(ax,[wavelength_nm(1),wavelength_nm(end)]);
ylim(ax,[0,1]);

end