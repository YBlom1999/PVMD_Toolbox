function S = Lux58(V,F,T)
%LUX determines the absorptivity in a user-defined object using ray-tracing
%and from this calculates the sensitivity S(absorptivity [W/m2] / beam
%normal irradiance [W/m2]).
%The object (e.g. a PV module) is built up from facets (triangles or
%parallelograms). The optical properties can be specified for each facet.
%---INPUT---
%V      vertex matrix. Each row has xyz coordinates of a vertex point.
%F      facet matrix. Each row has 3 or 4 vertex nrs specifying a facet.
%T      struct defining interface types (fields: Facet, RT, Emit, Plot)
%---OUTPUT---
%S      sensitivity
%b-version for building simulation (some redundant features disabled)
S = [];
%---check licence---
ok = license2kill(nargin);           %check licence
if ~ok || nargin == 0, return; end   %if licence not ok, don't continue

%---if 2 or less input arguments, generate preview and return---
if nargin==1, plot_V(V);   return; end     %if one input, assume V, plot vertices
if nargin==2, plot_F(V,F); return; end   %if two input, assume V and F

%---if all 3 input arguments are given, check if there are any rays to trace---
s = lux_settings2;                 %load settings
[s,T,Nr] = chkin(s,T);             %prepare parallel calculation (and accumulator)
if Nr == 0, plot_T(V,F,T); return; end %if no rays, just plot geometry with facets in true colors

%---when there are rays, LET THERE BE LIGHT!---
if s.plt==2 && ~s.parll, plot_T(V,F,T); end    %if rays are going to be plotted, plot geometry

TRI = reformat(V,F,T);           %calculate facet data matrix (normals etc)
[RAY,acz] = initial_ray(TRI,T,s);  %initialize rays

ACC = zeros(size(F,1)+3,s.nal+1,s.nw);  %accumulates energy absorbed in every
%absorber layer of every facet (+ weak/far/old)
%dim1:facet, dim2:layer, dim:wavelength

if s.pp                           %if parallel pool (distributed computing)
    parfor nr = 1:Nr              %for every ray
        acc = one_ray(RAY(nr,:),TRI,T,s); %fill 1-ray accumulator
        ACC = ACC + acc;            %add to total accumulator
    end
    %NOTE: by using cell or struct one could make separate accumulators for
    %for every surface type (with different dimensions). However, parfor
    %does not accept cell or struct, so need to keep a single accumulator.
else
    for nr = 1:Nr                 %for every ray
        acc = one_ray(RAY(nr,:),TRI,T,s); %fill 1-ray accumulator
        ACC = ACC + acc;           %add to total accumulator
    end
end
%NOTE: parfor for a large number of rays does not seem to create much
%overhead since the simulation time is reduced by a factor equal to the nr
%of cores. No need to optimize.

chkout(ACC,Nr)         %check output and generate warning if something suspicious

S = sensitive(ACC,TRI,Nr,acz,s);   %calculate sensitivity from accumulated irradiance

if s.plt==1, plot_S(V,F,T,S); end   %plot sensitivity of every facet in geometry

disp('')


%==========================================================================
    function ok = license2kill(n)
        %Compare current date to end of license date. If not yet expired ok = 1,
        %else ok = 0. Display message 2 weeks before expiration.
        
        licencee = 'Rudi Santbergen (TU Delft)';
        end_date = '25-May-2026';
        
        remain = datenum(end_date) - now;      %days of licence remaining
        
        if remain > 14                         %if licence valid for > 14 days
            ok = 1;                            %ok, no need to display message
            if n==0                            %if GenPro4 is run without input give general message
                disp('LUX ray-tracing model, by Rudi Santbergen, Delft University of Technology')
                disp(['Copy of ',licencee])
                disp(['Licence valid until: ',end_date])
            end
        elseif remain >= 0                     %if license valid for 0 to 14 days
            ok = 1;                            %ok, but display reminder
            disp('LUX ray-tracing model, by Rudi Santbergen, Delft University of Technology')
            disp(['Copy of ',licencee])
            disp(['Reminder: Your license will expire on ',end_date])
            disp('To extend your license please contact Rudi Santbergen (r.santbergen@tudelft.nl)');
        else
            ok = 0;                            %licence expired, display message
            disp('LUX ray-tracing model, by Rudi Santbergen, Delft University of Technology')
            disp(['Copy of ',licencee])
            disp(['Your license expired on ',end_date])
            disp('To renew your license please contact Rudi Santbergen (r.santbergen@tudelft.nl)');
        end
        
    end
%--------------------------------------------------------------------------
    function plot_V(V)
        %PLOT_V plots the vertex points in matrix V and labels them with
        %the corresponding vertex numbers.
        
        figure(1)
        clf
        %---plot vertices---
        line(V(:,1),V(:,2),V(:,3),'Marker','.','Color',[0.7 0.7 0.7],...
            'MarkerEdgeColor','r','MarkerSize',16);
        
        axis equal tight        %make xyz aspect ratio equal
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        zlabel('z (\mum)')
        view(30,30)            %standard view from the side
    end
%--------------------------------------------------------------------------
    function plot_F(V,F)
        %PLOT_F plots the facets in matrix F.
        
        figure(2)
        clf
        %---plot facets---
        p = patch('Faces',F,'Vertices',V);      %Matlab patch function
        set(p,'FaceAlpha',0.1,'FaceColor','r')  %make all transparent red
        
        axis equal off tight                    %xyz aspect ratio equal
        view(30,30)                            %standard view from side
    end
%--------------------------------------------------------------------------
    function plot_T(V,F,T)
        %PLOT_T plots the facets in matrix F (like plot_F) but gives each
        %facet the color specified by Type.
        
        figure(3)
        clf
        
        it = find(~cellfun(@isempty,{T.Invis})); %invisible (TODO: checks for not empty, value isn't considered)
        
        nrf = size(F,1);                         %nr of facets
        ix = true(nrf,1);
        for t = 1:length(it), ix(T(it(t)).Facet) = false; end  %set visibility 0 for invisible facets
        
        %---plot facets in their specified color---
        p = patch('Faces',F(ix,:),'Vertices',V);       %plot facets
        
        rgb = zeros(nrf,3);                      %initialize
        opa = zeros(nrf,1);                      %initialize
        for t = 1:length(T)
            if ~isempty(T(t).Facet)
                if isstruct(T(t).RT)
                    rgb(T(t).Facet,1) = 1;      %give cell red color
                    %TODO: somehow doesn't work
                else
                    rgb(T(t).Facet,:) = T(t).RT(1)*ones(numel(T(t).Facet),3);   %others gray-scale
                    opa(T(t).Facet) = 1-T(t).RT(2);
                end
            end
        end
        
        %set rgb and opacity
        set(p,'FaceVertexCData',rgb(ix,:),'CDataMapping','scaled',...
            'FaceColor','flat','FaceVertexAlphaData',opa(ix),...
            'AlphaDataMapping','none','FaceAlpha','flat')
        
        axis equal off tight                    %xyz aspect ratio equal
        view(30,30)                            %standard view from side
        compass(V)
        drawnow
    end
%--------------------------------------------------------------------------
    function compass(V)
        %display North-South-East-West in plot
        text(max(V(:,1)),               (min(V(:,2))+max(V(:,2)))/2,min(V(:,3)),'EAST');
        text((min(V(:,1))+max(V(:,1)))/2,max(V(:,2)),               min(V(:,3)),'NORTH');
        text(min(V(:,1)),               (min(V(:,2))+max(V(:,2)))/2,min(V(:,3)),'WEST');
        text((min(V(:,1))+max(V(:,1)))/2,min(V(:,2)),               min(V(:,3)),'SOUTH');
    end
%--------------------------------------------------------------------------
    function [s,T,Nr] = chkin(s,T)
        %CHKIN checks input argument T for consistency and start parpool if
        %needed.
               
        et = find(~cellfun(@isempty,{T.Emit}),1); %emitting type (max one)
        if isempty(et)                       %if no emitting type is defined
            Nr = 0;                          %there are no rays to trace
        else
            Nr = T(et).Emit(1);              %get number of rays
            if ~isfield(T,'Scat'), T(1).Scat = [];end       %make sure Scat field exists
            if ~isfield(T,'Teleport'), T(1).Teleport = [];end %make sure Teleport field exists
        end
        
        %---start parpool for parallel computation---
        willing = s.parll;                                       %is parallel pool desired? see lux_settings2.m
        able = license('test','Distrib_Computing_Toolbox');      %is Toolbox present?
        off = isempty(gcp('nocreate'));                          %is parallel pool turned off?
        
        if willing && able                                       %willing & able & off --> on
            s.pp = 1;
            if off, parpool; end
        else
            s.pp = 0;
        end
        %Note: Parpool is never turned off by LUX. But Matlab turns if off
        %automatically.
        
        %---check the maximum nr of absorber layers to size the accumulator---
        s.nal = 0;                           %assume 0 absorber layers
        s.nw = 1;                            %assume 1 wavelength
        for i = 1:length(T)                  %for each facet type
            if isstruct(T(i).RT)             %if is structure (CELL)
                sz = length(T(i).RT.lay);  %look at layer data
                s.nw = length(T(i).RT.wav);  %look at wavelength data
                %TODO:will cause problems if different facet types have
                %different wavelength vectors.
                
                %angle interpolation function (interpx1x) does not do
                %extrapollation. Need to set first 0deg and final 90deg,
                %which is not already the case if data comes from GP4.
                if T(i).RT.aoi(1)>0             %if first angle larger than 0
                    T(i).RT.aoi = [0;T(i).RT.aoi];            %add 0 to the front
                    T(i).RT.RAT = cat(2,T(i).RT.RAT(:,1,:),T(i).RT.RAT); %simply nearest neighbor
                end
                if T(i).RT.aoi(end)<90              %if final angle is less than 90
                    T(i).RT.aoi = [T(i).RT.aoi;90];            %add 90 to the end
                    T(i).RT.RAT = cat(2,T(i).RT.RAT,T(i).RT.RAT(:,end,:)); %simply nearest neighbor
                end
            else                             %if is matrix
                sz = size(T(i).RT,2);        %check size of RT matrix, nr of layers is evident from dim2.
            end
            s.nal = max(sz-2,s.nal); %update value when larger
            %Has at least 2 numbers: R and T.
            %Additional numbers correspond to absorber layers.
        end
    end
%--------------------------------------------------------------------------
    function TRI = reformat(V,F,T)
        %REFORMAT converts the input facet data (consisting of 3 or 4
        %vertex xyz coordinates) into a format suitable for intersection
        %test. Every row of TRI represents one facet and contains 19 nrs.
        % # 1: first vertex x-coordinate
        % # 2: first vertex y-coordinate
        % # 3: first vertex z-coordinate
        % # 4: edge vector E1 x-component
        % # 5: edge vector E1 y-component
        % # 6: edge vector E1 z-component
        % # 7: edge vector E2 x-component
        % # 8: edge vector E2 y-component
        % # 9: edge vector E2 z-component
        % #10: surface area of facet
        % #11: surface normal vector x-component
        % #12: surface normal vector y-component
        % #13: surface normal vector z-component
        % #14: a = E1*E1 (dot product pre-calculated)
        % #15: b = E2*E2 (dot product pre-calculated)
        % #16: c = E1*E2 (dot product pre-calculated)
        % #17: d = c^2-ab (denominator pre-calculated)
        % #18: is triangle? (1=triangle (3angle), 0=parallelogram (4angle))
        % #19: type number (for finding correct optical properties)
        
        %---both 3 and 4angle can be characterized by 3 vertices in FF ---
        nf = size(F,1);                       %nr of facets
        FF = zeros(nf,3);                     %initialize FF
        is3 = zeros(nf,1);                    %distinguishes 3 and 4 angle
        for f = 1:nf                          %for every row of F
            ix = ~isnan(F(f,:));              %indicate non-NaNs
            if sum(ix)==4                     %if 4angle
                FF(f,:) = F(f,[1 2 4]);       %take first, second and final
            else                              %if 3angle
                FF(f,:) = F(f,ix);            %take all non-NaNs
                is3(f) = 1;                   %remember that it's 3angle
            end
        end
        
        %---FF was created because below is convenient---
        V0 = V(FF(:,1),:);          %first vertex
        V1 = V(FF(:,2),:);          %first vertex' one side neighbor
        V2 = V(FF(:,3),:);          %first vertex' other side neighbor
        
        E1 = V1-V0;                 %edge vector1
        E2 = V2-V0;                 %edge vector2
        
        %---pre-calculation of dot products---
        a = E1(:,1).*E1(:,1)+E1(:,2).*E1(:,2)+E1(:,3).*E1(:,3);  %E1*E1
        b = E2(:,1).*E2(:,1)+E2(:,2).*E2(:,2)+E2(:,3).*E2(:,3);  %E2*E2
        c = E1(:,1).*E2(:,1)+E1(:,2).*E2(:,2)+E1(:,3).*E2(:,3);  %E3*E3
        
        d = c.*c-a.*b;      %pre-calculation of denominator
        
        %---edge vector cross-product gives surface normal vector---
        Nx = E1(:,2).*E2(:,3)-E1(:,3).*E2(:,2);     %x-component
        Ny = E1(:,3).*E2(:,1)-E1(:,1).*E2(:,3);     %y-component
        Nz = E1(:,1).*E2(:,2)-E1(:,2).*E2(:,1);     %z-component
        %sign convention such that normal points to RIGHT HAND RULE side
        
        Ap = sqrt(Nx.^2+Ny.^2+Nz.^2);     %geometrical length of N <---This is the surface area of parallelogram
        Nx = Nx./Ap;                      %normalize xyz components
        Ny = Ny./Ap;                      %(important for calculation
        Nz = Nz./Ap;                      %of direction R&T rays.)
        
        A = Ap;                            %facet surface area is parallelogram area
        A(logical(is3)) = Ap(logical(is3))/2;  %except triangle has half the surface area
        
        tv = zeros(nf,1);                %initialize type vector
        for t = 1:length(T)              %for every type
            ix = T(t).Facet;             %corresponding facets nr
            tv(ix) = t;                  %vector element will be type
        end
        
        %---combine everything in TRI matrix---
        TRI = [V0,E1,E2,A,Nx,Ny,Nz,a,b,c,d,is3,tv];
        
    end
%--------------------------------------------------------------------------
    function [RAY,acz] = initial_ray(TRI,T,s)
        %INITIAL_RAY creates the initial ray matrix of which each row
        %represents one ray emitted by the light source. Each ray is
        %characterized by 9 numbers (i.e. matrix has 9 columns)
        % #1: ray origin x-coordinate
        % #2: ray origin y-coordinate
        % #3: ray origin z-coordinate
        % #4: ray direction x-component
        % #5: ray direction y-component
        % #6: ray direction z-component
        % #7: ray generation nr. (nr. bounces)
        % #8: ray plot (0 = no, 1 = yes)
        % #9: ray intensity [-]
        
        et = find(~cellfun(@isempty,{T.Emit}),1); %emitting type (just one)
        ef = T(et).Facet;             %emitting facets (can be multiple)
        a = sum(TRI(ef,10));          %total surace area of emitting facet(s)
        emit = T(et).Emit;            %get input from struct
        
        [ray_O,Nrr] = ray_origin(emit(1),TRI,ef);      %ray origin
        
        [ray_D,theta] = ray_direction(emit(2:3),Nrr);          %ray direction
        acz = a*cosd(theta);       %projected area of light source (for calculation of sensitivity)
        %NOTE: it is implicitely assumed that the light source is
        %horizontal, radiating downward!!!!!!!!
        %TODO: maybe generate warning when this is not the case
        
        NrWav = s.nw;                   %set nr of wavelengths
        ray_I = ones(Nrr,NrWav);        %ray intensity
        
        %---ray generation (nr. interfaces touched)---
        ray_G = zeros(Nrr,1);           %all rays start as generation 0
        
        %---ray plotability (plot takes time, not all rays are plotable)---
        ray_plt = zeros(Nrr,1);         %initialize vector
        nr_plot = ceil(sqrt(Nrr));      %sqrt rule (100 rays, 10 plotable)
        ix = randperm(Nrr,nr_plot);     %select plotable rays randomly
        ray_plt(ix) = 1;                %make selected rays plotable
        
        %---combine all information in one matrix---
        RAY = [ray_O,ray_D,ray_G,ray_plt,ray_I];   %nr. rows = nr. of rays
        %..................................................................
        function [ray_O,Nrr] = ray_origin(emit1,TRI,ef)
            %determine ray origin xyz coordinates based on input1 (first
            %cell of emit input-field)
            
            Nrr = emit1;                     %represents number of rays
            nrr = ceil(Nrr/length(ef));      %nr of rays per facet
            Nrr = length(ef)*nrr;            %total nr of rays (in reality)
            of = ones(nrr,1);                %ones to get dimensions right
            o3 = ones(1,3);                  %ones to get dimensions right
            
            ray_O = zeros(Nrr,3);            %initialize ray origins xyz
            fix = 1:nrr;                     %facet index
            for f = ef                       %for every emitting facet
                [r1,r2] = quasirandom(nrr);  %quasi random numbers
                
                X1 = of*TRI(f,1:3);          %first vertex point
                E1 = of*TRI(f,4:6);          %edge vector 1
                E2 = of*TRI(f,7:9);          %edge vector 2
                O1 = X1+(r1*o3).*E1+(r2*o3).*E2;        %points on 4angle
                O2 = X1+((1-r1)*o3).*E1+((1-r2)*o3).*E2;%same points mirrored
                ixa = 1:nrr;                 %index of all points
                if TRI(f,18), ixa  = r1+r2<1; end %logical index inside points
                %outside 3angle points are replaced by their mirror point
                ray_O(fix,:) = [O1(ixa,:);O2(~ixa,:)];
                fix = fix + nrr;             %increment for next emitting facet
            end
            
            %............................................................
            function [X,Y] = quasirandom(nr)
                %generates set of nr 2D coordinates by creating a regular grid
                %with a point at a random position in each grid cell. This
                %gives a more uniform distribution compared to  truly random
                %points. Note: Sobolset could be used but requires statistics
                %toolbox
                
                q = ceil(sqrt(nr));     %nr of subcells in 1D
                q2 = q^2;               %nr of subcells in 2D (has excess)
                cw = 1/q;               %subcell width
                x = 0:cw:(1-cw);        %subcell centre grid coord (x=y)
                
                %---similar to meshgrid but creates vector (not matrix)---
                Xg = zeros(q2,1);       %initialize grid coordinate x
                Yg = zeros(q2,1);       %initialize grid coordinate y
                
                for cc = 1:q            %for every row
                    ixx = (cc-1)*q+(1:q); %row index
                    Xg(ixx) = x;        %all grid coordinates (x varies)
                    Yg(ixx) = x(cc);    %one grid coordinate (y is constant)
                end
                %---
                
                Xr = cw*rand(q2,1);     %random component x
                Yr = cw*rand(q2,1);     %random component y
                
                X = Xg+Xr;              %combine grid and random components x
                Y = Yg+Yr;              %combine grid and random components y
                
                ex = q2-nr;             %q2-nr is excess points
                for cc = 1:ex           %for every excess point
                    d = randi(length(X),1); %select random element
                    X(d) = [];          %remove it from X
                    Y(d) = [];          %remove it from Y
                end
            end
        end
        %..................................................................
        function [ray_D,theta] = ray_direction(emit23,Nrr)
            %Determine ray direction xyz coordinate based on second and
            %third input in Emit vector.
            %Azimuth angle: 0 = N, 90 = E, 180 = S, 270 = W.
            
            theta = emit23(1);   %zenith angle [deg]
            Theta = theta*ones(Nrr,1);
            phi = emit23(2);     %azimuth [deg]
            Phi = phi*ones(Nrr,1);
            
            rho = sind(theta);
            Dx = rho.*sind(Phi);            %direction vector x component
            Dy = rho.*cosd(Phi);            %direction vector y component
            Dz = -cosd(Theta);              %direction vector z component
            ray_D = [Dx,Dy,Dz];             %ray direction xyz-list
            
            %Note: zenith angle theta is defined relative to -z axis, not
            %relative to source normal!
            %Note: azimuth 0 deg is south. This means rays are coming FROM 
            %south but are pointing TO north.
        end
        %..................................................................
    end
%--------------------------------------------------------------------------
    function chkout(ACC,Nr)
        %check output and generate warning if something suspicious
        
        %Every ray has 1 unit of energy. All energy should be accounted for.
        Ea = squeeze(sum(ACC(:,1,:),1));   %energy absorbed, total in all facets (wavelength vector)
        dE = max(abs((Nr-Ea)/Nr));         %relative energy loss (max value)
        if dE>1e-6, warning([num2str(dE),' of energy unaccounted for']); end
        
        %last three entries in ACC do not correspond to a facet
        far = ACC(end-2,1,:)/Nr;
        old = ACC(end-1,1,:)/Nr;
        weak = ACC(end,1,:)/Nr;
        
        if max(far)>1e-6, warning('Rays have escaped the simulation domain'); end
        if max(old)>1e-3, warning('Old rays terminated. Increase ''S.mx_g''.'); end
        if max(weak)>5e-3, warning('Weak rays terminated. Decrease ''S.Imin''.'); end
        
    end
%--------------------------------------------------------------------------
    function S = sensitive(ACC,TRI,Nr,acz,s)
        %calculate sensitivity from accumulated irradiance
        %INPUT
        %ACC        accumulated absorbed intensity [W] (dim1: facet, dim2: layer (tot,A1,A2,...), dim3: wav)
        %TRI        facet information
        %Nr         number of rays. Each ray has unit intensity (1W)
        %acz        cosine for projected area of light source
        %s          settings
        %OUTPUT
        %S          sensitivity with respect to absorbed power [-] (dim1: facet, dim2: layer)
        ACC(end-2:end,:,:) = [];    %remove far/weak/old accumulators (not linked to any facet)
        A = TRI(:,10);              %facet area
        %S = (acz/Nr)*ACC./(A*ones(1,s.nal+1)); %sensitivity with respect to ABSORBED irradiance
        S = (acz/Nr)*ACC./repmat(A,1,s.nal+1,s.nw); %sensitivity with respect to ABSORBED irradiance
    end
%--------------------------------------------------------------------------
    function plot_S(V,F,T,S)
        %plot sensitivity of every facet in geometry
        
        ix = true(size(F,1),1);                  %set visible is true for all facets
        it = find(~cellfun(@isempty,{T.Invis})); %find invisible types
        for t = 1:length(it), ix(T(it(t)).Facet) = false; end %set visible to false for all invisitble facets
        
        et = find(~cellfun(@isempty,{T.Emit}),1);
        
        %for l = 1:size(S,2)     %for every layer in multilayer structure
        for l = 1:1              %plot only the 'total' absorptance sensitivity
            %TODO: will show zeros for single layer surfaces when plotting
            %second layer.
            figure(4+l)
            clf
            patch('Faces',F(ix,:),'Vertices',V,'FaceVertexCData',S(ix,l),'FaceColor','flat')
            %includes only the visible facets
            caxis([0,1])
            colorbar
            shading flat
            axis equal off tight                    %xyz aspect ratio equal
            compass(V)
            
            nrr = T(et).Emit(1);                %nr of rays
            zenith  = round(T(et).Emit(2));     %zenith angle
            azimuth = round(T(et).Emit(3));     %azimuth angle
            title([num2str(nrr),' rays',10,'zenith: ',num2str(zenith),'\circ',10,'azimuth: ',num2str(azimuth),'\circ'])
            
            view(-azimuth,90-zenith)                            %standard view from side
            drawnow
            %pause
        end
    end
%--------------------------------------------------------------------------
end

%Functions below should not be nested to allow parallel processing
%==========================================================================
function [acc] = one_ray(ray,TRI,Type,s)
%Traces one ray (including descendants) until all its energy is
%accumulated
%INPUT
%ray        incoming ray origin, direction and intensity (TODO: wavelength
            %dependend intensity--> extend vector) 
%TRI        facet information
%Type       facet types
%s          ray-tracing settings
%OUTPUT
%acc        accumulator (dim1: facet nr., dim2: layer nr., dim3:
%           wavelength)

acc = zeros(size(TRI,1)+3,s.nal+1,s.nw);    %initialize intensity accumulator
%dim1 is equal to the nr of facets + 3 (old/weak/far)
%dim2 is equal to the nr of absorber layers + 1. The first column
%       contains the segment's total absorptance (1-R-T) which is needed for
%       energy conservation check
%dim3 is equal to the nr of wavelengths (TODO: for now set to 1)

RAYQ = zeros(10,size(ray,2));          %initialize ray-queue (length 10)
%Note:  Initial queue length of 10 is more than sufficient. Only when rays are
%wildly branching queue is expanded, which might take some time
%TODO: There is still a problem with index going to 0 when transmitted rays
%are involved (transparent surfaces).

RAYQ(1,:) = ray;               %add the original ray to ray-queue
active = 1;                    %pointer to the column of the active ray
while active>=1                %while there is an active ray
    ray_in = RAYQ(active,:);   %take active ray from queue
    [ray_out,acc] = one_segment(ray_in,TRI,Type,s,acc); %let it bounce
    nr_out = size(ray_out,1);  %nr of outgoing rays
    if nr_out>0, RAYQ(active:active-1+nr_out,:)=ray_out; end %add outgoing rays to queue (if any, none when ray hits perfect absorber)
    active = active-1+nr_out;  %increment pointer accordingly (do as last)
    %such that next ray on queue becomes active (one ray done, nr_out rays
    %added)
end

%---Energy conservation check---
% dE1 = max(abs(squeeze(sum(acc(:,1,:),1))-1));       %test energy conservation
% if dE1>1e-6, warning([num2str(dE1),' of energy lost of when tracing single ray']);end

%--------------------------------------------------------------------------
    function [ray_out,acc] = one_segment(ray_in,TRI,Type,s,acc)
        
        %---for energy conservation test---
        %E00=sum(ray_in(:,9:end),1)';
        %E01=squeeze(sum(acc(:,1,:)));      %energy content of accumulator at start wavelenght vector
                
        I_in = ray_in(9:end);              %incident ray intensity
        %---calculate ray's intersection point---
        [X,N,i] = intersexy(ray_in,TRI);   %get intersection point X, local normal N and facet index i
        t = TRI(i,19);                  %intersected facet type
               
        if isempty(X), acc(end-2,1,:) = acc(end-2,1,:)+shiftdim(I_in,-1); ray_out =[]; return; end %far ray accumulation
        
        if s.plt==2 && ray_in(8), figure(3); line([ray_in(1),X(1)],[ray_in(2),X(2)],[ray_in(3),X(3)],'Color','w'); drawnow; end %plot segment (origin to intersection)
        
        if ray_in(7)>s.mx_g, acc(end-1,1,:) = acc(end-1,1,:)+shiftdim(I_in,-1);ray_out = []; return; end   %old ray accumulation
        
        [ray_R,ray_T,I_a] = reflexy(ray_in,X,N,Type(t),s);
        
        %---facet ACCUMULATION!!!---
        acc(i,1:size(I_a,2),1:size(I_a,1)) = acc(i,1:size(I_a,2),1:size(I_a,1))+shiftdim(I_a',-1);                   
        %The accumulator is a matrix, containing the ABSORBED irradiance
        %for every facet. Some facets represent multilayer systems and then
        %the absorptance in every absorber layer is accumulated. The number
        %of absorber layers can vary for different facet types. In that
        %case the left-most columns are filled.
        
        if mean(ray_R(9:end))<s.Imin, acc(end,1,:) = acc(end,1,:)+shiftdim(ray_R(9:end),-1); ray_R = []; end   %weak ray accumulation
        if mean(ray_T(9:end))<s.Imin, acc(end,1,:) = acc(end,1,:)+shiftdim(ray_T(9:end),-1); ray_T = []; end   %weak ray accumulation
        %Note: ray is weak when its wavelength AVERAGE intensity falls
        %below the threshold value.
        
        if ~isempty(Type(t).Scat)
            [ray_R] = surf_scat(ray_R,N,Type(t).Scat,s);     %surface scattering (reflected ray only)
        end
        
        ray_out = [ray_R;ray_T];
        
        %---test energy conservation for debugging---
        %E02 = squeeze(sum(acc(:,1,:)));   %energy content of accumulator at end (wavelength vector)
        %E03 = sum(ray_out(:,9:end),1)';
        %if numel(E03)==0, E03 = 0; end
        
        %dE0 = E00+E01-E02-E03;  %energy in minus out     
        %if max(abs(dE0))>1e-6
        %    warning([num2str(max(abs(dE0))),' of energy lost when tracing single ray-segment']);
        %end
        
        %..................................................................
        function [X,N,i] = intersexy(ray,TRI)
            %Calculates nearest intersection between 1 ray and all
            %triangles.
            %INPUT
            %ray        single ray information (see initial_ray function)
            %TRI        all triangle information (see reformat function)
            %OUTPUT
            %X          intersection coordinate xyz
            %i          index of intersected triangle
            
            nr3 = size(TRI,1);              %nr. of triangles
            %---ray data (from initial_ray function)---
            P0 = ones(nr3,1)*ray(1:3);
            D  = ones(nr3,1)*ray(4:6);
            %---triangle data (from reformat function)---
            V0 = TRI(:,1:3);
            E1 = TRI(:,4:6);
            E2 = TRI(:,7:9);
            NN = TRI(:,11:13);
            a  = TRI(:,14);
            b  = TRI(:,15);
            c  = TRI(:,16);
            d  = TRI(:,17);
            is3= TRI(:,18);
            
            q = dotx(NN,(V0-P0))./dotx(NN,D); %intersec. distance (units D)
            
            Pi = P0+(q*ones(1,3)).*D;         %potential intersection point for every triangle
            
            E3 = Pi-V0;                       %vertex-to-intersection_point vector
            dE13 = dotx(E1,E3);               %dot product E1*E3
            dE23 = dotx(E2,E3);               %dot product E2*E3
            
            rr = (c.*dE23-b.*dE13)./d;         %units of E1 to reach intersection point
            ss = (c.*dE13-a.*dE23)./d;         %units of E1 to reach intersection point
            
            %---intersection test---
            ix = (q>1e-9 & rr>=0 & ss>=0 & ((is3 & rr+ss<=1) | (~is3 & rr<=1 & ss<=1)));
            %note q>1e-9 to avoid self-intersection
            %TODO: creates escape possibility (theoretically, no problems
            %in practice), alternative test?
            
            if sum(ix)==0                     %if none pass (shouldn't happen)
                i = [];                       %no intersected triangle
                X = [];                       %no intersection point
                N = [];                       %no surface normal
            else                              %if one or more pass test
                q(~ix) = Inf;                 %assign Infinity to not passed test
                [~,i] = min(q);               %find closest among passed
                X = Pi(i,:);                  %corresponding intersection xyz
                N = NN(i,:);                  %corresponding surface normal xyz
            end
            
        end
        %..................................................................
        function [ray_R,ray_T,I_a] = reflexy(ray_in,X,N,Type,s)
            %function [ray_R,ray_T,ACi1,ACr1,ACf1,pr1] = reflexy(ray_in,ACi1,ACr1,ACf1,pr1,X,N,Tw,TRI,fx,i,S)
            %calculates reflected and transmitted ray from incident ray.
            %Absorbed energy is accumulated.
            
          
            Di = ray_in(4:6);           %direction incident ray
            cosphi = Di*N';             %cos(angle of inc)
            
            %---properties of reflected and transmitted rays---
            Dr = Di-2*cosphi*N;         %direction reflected ray
            
            [Ot] = origint(Type,X,N);  %origin transmitted ray (it may teleport)
            
            %---ANGULAR DEPENDENCE IS INCLUDED HERE (new in this version)---
            %could convert lookup table from phi to cosphi to avoid acosd
            %inside the loop. But not worth it as acosd hardly takes up
            %time.
            
            if isstruct(Type.RT)  %RT is a struct containing 3D array for wavelength, angle and layer dependence
                %angle dependent, so need to get phi
                phi = acosd(cosphi);                                      %is between 0 and 180 deg!
                if phi > 90, phi = 180-phi; end                           %keep phi between 0 and 90 deg
                
                %Type.RT.RAT = 46 x 15 x 4
                %after interpolating the angle
                %RAAT = 46 x 4
                %first column spectral reflectance 46 x 1
                %final column spectral transmittance 46 x 1
                %middle columns spectral absorptance 46 x 3 (total and per absorber layer!)
                
                RAAT = interpx1x(Type.RT.aoi,Type.RT.RAT,phi);  %use 3D interpolator x1x
                R = RAAT(:,1);
                T = RAAT(:,end);
                A = [1-R-T,RAAT(:,2:end-1)]; 
            else                     %if RT is a matrix
                sz = size(Type.RT);
                o = ones(s.nw,1);
                if sz(1)==1          %if RT is a row vector
                    R = Type.RT(1)*o;        %take fixed values(angle and wavelength independent)
                    A = [1-Type.RT(1)-Type.RT(end),Type.RT(2:end-1)]*o;
                    T = Type.RT(end)*o;
                else                 %if it really is a 2D matrix
                    %angle dependent, so need to get phi
                    phi = acosd(cosphi);                                      %is between 0 and 180 deg!
                    if phi > 90, phi = 180-phi; end                           %keep phi between 0 and 90 deg
                    
                    RAAT = interp1x(Type.RT(:,1),Type.RT(:,2:end),phi);
                    R = RAAT(1)*o;            %reflectance
                    T = RAAT(end)*o;          %transmittance
                    A = [1-R-T,RAAT(2:end-1)]*o; %total absorptance followed by layer absorptances
                end
            end
            %TODO: does this all have to be done on the fly, or can some
            %precalculation be done for more speed?
            %---
            %total absorptance
            
            I_i = ray_in(9:end)';  %intensity incident ray (function of wavelength)
            I_r = I_i.*R;          %intensity reflected ray (function of wavelength)
            I_t = I_i.*T;          %intensity transmitted ray (function of wavelength)
            o = ones(1,size(A,2));
            I_a = (I_i*o).*A;      %intensity absorbed at facet (total, lay1, lay2,...)
                                  %function of both layer and wavelength
            
            ray_R = [X ,Dr,ray_in(7)+1,ray_in(8),I_r']; %reflected ray (gen+1,plot same)
            ray_T = [Ot,Di,ray_in(7)+1,ray_in(8),I_t']; %transmitted ray
            %TODO: should be possible with less transpositions
            
            dE = I_i - I_r - I_t - I_a(:,1);        %energy in minus out
            if max(abs(dE))>1e-6, warning([num2str(max(abs(dE))),' of energy lost of in single reflection']); end
            
            %..............................................................
            function Ot = origint(Type,X,N)
                %Determine origin of transmitted ray. Normally it is the
                %intersection point, but in case of teleport (or
                %randomizing plane) it is not.
                
                if isempty(Type.Teleport)   %if not defined
                    Ot = X;               %origin is intersection point (trivial)
                else
                    Ot = X + Type.Teleport*N; %teleport in normal direction
                end
            end
            %..............................................................
            function yi=interp1x(xjo,yjo,xi)
                %super-mega fast linear interpolation (modified interp1)
                %interpolates along the first dimension of a 2D array 
                
                r=find(xjo<=xi,1,'last');           %index of previous (next is previous + 1)
                u = (xi-xjo(r))/(xjo(r+1)-xjo(r));  %relative position of xi between previous and next
                yi=yjo(r,:)+(yjo(r+1,:)-yjo(r,:))*u;
            end
            %..............................................................
            function yi=interpx1x(xjo,yjo,xi)
                %super-mega fast linear interpolation (modified interp1)
                %interpolates along the second dimension of a 3D array
               
                r=find(xjo<=xi,1,'last');           %index of previous (next is previous + 1)
                u = (xi-xjo(r))/(xjo(r+1)-xjo(r));  %relative position of xi between previous and next
                yi=yjo(:,r,:)+(yjo(:,r+1,:)-yjo(:,r,:))*u;
            end
            
        end
        %..................................................................
        function [ray_out] = surf_scat(ray_in,N,Scat,S)
            %Scatter one incident ray into several outgoing rays. Note:
            %scattering is applied AFTER reflection. ray_in is reflected
            %ray
            %TODO: this function is rather time consuming. Much could be
            %done outside the loop.
            
            if isempty(ray_in), ray_out = []; return; end
            
            H = Scat(1);                  %haze value
            ph = Scat(2);                 %AID parameter (for now: Phong exponenent)
            
            %---specular ray----
            if numel(H)>1 || H<0 || H>1
                error('Unrecognized Haze value');
            else
                if H == 0
                    ray_out = ray_in;                   %pure specular reflection (better to remove scat-field from input)
                    return
                elseif H==1                             %if haze == 1,
                    ray_spec = [];                      %no specular ray
                else
                    ray_spec = ray_in;                  %same as incommning ray
                    ray_spec(9) = ray_spec(9)*(1-H);    %intensity reduced by factor (1-H)
                end
            end
            
            %---diffuse rays---
            %TODO: one could pre-calculate a set of diffuse rays, store
            %them in Scat and use them here. This would save some
            %computation time.
            
            I_dif = H*ray_in(9:end);                        %total intensity of diffuse rays
            n_dif = ceil(mean(I_dif)/S.Imc);                  %nr of diffuse rays.
            %below Imc -->1 diffuse ray
            %above Imc -->multiple diffuse rays with intensity on the order of mc-threshold
            if n_dif>0
                [ang1,ang2] = randang(n_dif,ph);        %inclination and azimuth angles
                %base ray will be tilted by inclination angle and rotated
                %by azimuth angle.
                sg = sign(dotx(ray_in(4:6),N));         %side of incidence +1 = normal, -1 = opposite (ray_in is reflected ray! Parallel to N?)
                ray_dif = ones(n_dif,1)*ray_in;         %copy vector to matrix
                for a = 1:n_dif
                    V = rotate(sg*N,ang1(a),ang2(a)); %rotated direction vector
                    %V = rotate(ray_in(4:6),ang1(a),ang2(a)); %rotated direction vector
                    
                    %Note: 1st option distributes diffuse rays around N
                    %(minus N when reflected off rear side )
                    %2nd option around specular reflected ray (true PHONG
                    %style). This causes complications for rays rotated beyond
                    %interface. Therefore 2nd option is disabled for
                    %now. Maybe no problem since there still is the
                    %option to use the specular beam by setting H<1.)
                    
                    ray_dif(a,4:6) = V;
                    ray_dif(a,9:end) = I_dif/n_dif;     %ray intensity
                end
            else
                ray_dif = [];
            end
            
            ray_out = monte_carlo([ray_spec;ray_dif],S);
            
            %..............................................................
            function [ang1,ang2] = randang(nr,ph)
                %set of random angles. ang1 0-90 range, according to
                %probability distribution (Phong distribution). ang2 0-360
                %range, distributed uniformly.
                ang1 = probably(2*nr,ph);                %inclination angle
                ang2 = 360*rand([2*nr,1]);               %azimuth angle
            end
            %..............................................................
            function x = probably(nr,ph)
                %generate random numbers in 0-90 range according to Phong
                %distribution.
                xr = 90*rand(100*nr,1);                  %random nrs 0-90
                yr = rand(100*nr,1);                     %random nrs 0-1
                y_acc = PD(xr,ph);                       %y-acceptance level for every xy pair
                ixy = find(yr<y_acc,nr,'first');         %index of accepted xy pairs
                x = xr(ixy);                             %accepted x
                
                function y = PD(x,ph)
                    %probability distribution (normalized 0-1)
                    y = cosd(x).^ph;                %Phong distr. (=cos^ph)
                end
            end
            %..............................................................
            function V3 = rotate(V1,ang1,ang2)
                %rotate vector V1 over inclination angle ang1 and azimuth
                %angle ang2
                
                J = rand(1,3);                  %random vector
                ax1 = cross(J,V1);              %inclination axis
                ax1 = ax1/sqrt(sum(ax1.^2));
                V1 = V1/sqrt(sum(V1.^2));
                
                M1 = rotmatrix(ax1,ang1*pi/180);
                V2 = M1*V1.';
                
                M2 = rotmatrix(V1,ang2*pi/180);
                V3 = M2*V2;
                %..........................................................
                function M = rotmatrix(ax,ang)
                    %calculate 3x3 rotation matrix for rotation around axis
                    %ax over angle ang
                    x = ax(1);y = ax(2);z = ax(3);
                    cs = cos(ang);sn = sin(ang);C = 1-cs;
                    
                    %----definition of rotation matrix in 3D---
                    M = [x*x*C+cs   , x*y*C-z*sn , x*z*C+y*sn;
                        y*x*C+z*sn  , y*y*C+cs   , y*z*C-x*sn;
                        z*x*C-y*sn  , z*y*C+x*sn , z*z*C+cs   ];
                end
            end
            %..............................................................
            function ray_out = monte_carlo(ray_out,S)
                %If the total intensity of the outgoing rays is below the mc
                %threshold, one outgoing ray is selected which gets all intensity,
                %all others are terminated. This prevents the number of rays from
                %exploding.
                
                if isempty(ray_out), return; end %nothing left to do
                Imc = S.Imc;                %monte-carlo threshold
                In = ray_out(:,7);          %individual ray intensities
                I = sum(In);                %total intensity of outgoing rays
                nr = size(ray_out,1);       %nr of outgoing rays
                if I < Imc && nr>1          %if total intensity low and multiple rays
                    rnd = I*rand(1);        %random nr between 0 and I
                    for c = 1:nr            %for every outgoing ray
                        rnd = rnd - In(c);   %reduce random nr by ray intensity
                        if rnd<=0                   %if it becomes negative
                            ray_out = ray_out(c,:); %select current ray, remove others
                            ray_out(7) = I;   %current ray gets all intensity
                            break                   %job done, break from inner for-loop
                        end
                    end
                end
            end
            %..............................................................
        end
        %..................................................................
        function C = dotx(A,B)
            %dot product along dim2 of matrices A and B
            %faster than standard dot product
            
            C = A(:,1).*B(:,1)+A(:,2).*B(:,2)+A(:,3).*B(:,3);
        end
        %..................................................................
    end
end
