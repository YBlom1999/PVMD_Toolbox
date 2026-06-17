function [Rpa,Rsa,Rpb,Rsb,Apa,Asa,Apb,Asb,Tpa,Tsa,Tpb,Tsb,AngleTree,logbook]= RAY35at(N1,N2,Int,SET)

%note: 75% of angle tree is empty. Use sparse matrix instead? No, this is
%only possible for 2D matrices. AT is 3D.
%TODO: Many rays follow similar path. Use clever compression method
%combining similar rays.
%TODO: air/Si and Si/air interface are calculated separately, but have
%'flipped' AT. Fliping possible, but a bit complex in combination with skip
%dark (could be alternative for skip dark, be gain is limited)

Pre=PRE(N1,N2,Int.coat,SET.wavelength_um(SET.w));   %construct precalculation table

%---part A: calculate the ray paths---
%if isfield(Int,'AT') && ~isempty(Int.AT{SET.w})  %if precalculated AT given as input
if isfield(Int,'AT') && ~isempty(Int.AT)  %if precalculated AT given as input
    %TODO: check this. Might go wrong if AT does not have proper format.
    AngleTree = Int.AT{SET.w};              %use that
else                                      %otherwise...
    AngleTree = RAY35_A(Pre,Int,SET);     %...calculate using ray-tracing           
    %note: this part does not change when the coating is modified, will be
    %skipped if AT already pre-calculated
end

%---part B:calculate ray intensities---
[Rpa,Rsa,Rpb,Rsb,Apa,Asa,Apb,Asb,Tpa,Tsa,Tpb,Tsb,logbook]= RAY35_B(AngleTree,Pre,SET);
%Fresnel equation is used only here

if ~SET.out_ang_tree, AngleTree = []; end    %if not needed don't store AngleTree

% PLT(Rpa,Rsa,Tpa,Tsa,Rpb,Rsb,Tpb,Tsb)
% title(S.wav(S.w))
% drawnow


%--------------------------------------------------------------------------
      function Pre=PRE(N1,N2,coat,w)
        %create precalculation table for ray-tracing
        
        %INPUT
        %N1,N2      refractive indices on either side of interface
        %coat       struct with coating info (refr. index, thickness)
        %OUTPUT
        %Pre        precalc. table (double sided, incl. multi coat abs)
        
        eps=1e-9;       %small nr.
        nang=90;        %nr. of precalc angles (reciprocity!!!)
        aiB=linspace(0,pi/2,nang).'; %angle of incidence (FOR TABLE ONLY!)
        aiB(end)=[];    %exclude grazing inc.
        cai=cos(aiB);   %cos(angle of inc.)
        cai(1)=1+eps;   %avoid lookup on the edge
        
        %---create input for 'fresnel' model---
        cnr = length(coat);     %nr of coatings
        scnr = 0;               %initialize nr of SUBcoatings
        NN=zeros(2+cnr,1);
        tt=cell(2+cnr,1); %preallocations
        NN(1)=N1;
        tt{1}=inf;
        for c=1:cnr
            NN(1+c)=coat(c).N(SET.w);
            tt{1+c}=coat(c).thi;  %<--- this can be thickness VECTOR
            scnr = scnr + max(1,length(coat(c).thi)-1);
        end
        NN(2+cnr)=N2;
        tt{2+cnr}=inf;
        
        Pre1=zeros(nang,7+2*scnr);Pre2=Pre1;  %initialize table (2 halves)
        %===1st half: N1-->N2===
        %---thin-film optics model to calculate interface R, A, T
        [RATp,RATs,arr]=fresnel06(NN,tt,aiB,w);
        
        Pre1(1:nang-1,:)=[cai,cos(arr),sin(arr),RATp',RATs']; %TABLE 1
        
        %---add special case: grazing incidence---
        if N2<=N1;crit=pi/2;else, crit=asin(N1/N2);end
        RATg=[1,zeros(1,1+scnr)]; %grazing --> R=1
        Pre1(nang,:)=[eps,cos(crit),sin(crit),RATg,RATg]; %final row
        %---correct special case: TIR---
        for c=1:nang
            if Pre1(c,2)<eps; Pre1(c,2:end)=[0,1,RATg,RATg];end
        end
        Pre1=real(Pre1);    %for safety
        
        %===2nd half: N2-->N1===
        [RATp,RATs,arr]=fresnel06(flipud(NN),flipud(tt),aiB,w);
        %---break up RAT into R, A and T and flip the A---
        Rp=RATp(1,:)';
        Ap=RATp(2:1+scnr,:)'; Ap=fliplr(Ap);
        Tp=RATp(end,:)';
        Rs=RATs(1,:)';
        As=RATs(2:1+scnr,:)'; As=fliplr(As);
        Ts=RATs(end,:)';
        %---
        Pre2(1:nang-1,:)=[-cai,cos(arr),sin(arr),Rp,Ap,Tp,Rs,As,Ts]; %TABLE 2
        
        %---add special case: grazing incidence---
        if N2>=N1; crit=pi/2; else, crit=asin(N2/N1); end
        Pre2(nang,:)=[-eps,cos(crit),sin(crit),RATg,RATg];
        %---correct specal case: TIR---
        for c=1:nang
            if Pre2(c,2)<eps; Pre2(c,2:end)=[0,1,RATg,RATg];end
        end
        Pre2=real(Pre2);
        %===
        Pre=[Pre1;flipud(Pre2)];  %COMBINE table 1 & 2
    end
%--------------------------------------------------------------------------

%     function PLT(Rpa,Rsa,Tpa,Tsa,Rpb,Rsb,Tpb,Tsb)
%         %plot geometery and S-matrices (p&st pol)
%         
%         %---plot S-matrix (p-pol)---
%         figure(2)
%         RT4p=[Rpa,fliplr(Tpa);flipud(Tpb),rot90(Rpb,2)];    %join matrices
%         RT4s=[Rsa,fliplr(Tsa);flipud(Tsb),rot90(Rsb,2)];    %join matrices
%         imagesc(RT4p+RT4s)                                       %plot
%         axis equal tight xy
%         colormap hot
%         caxis([0,0.2])
%         set(2,'Position',[800 200 500 400])
%         drawnow
%     end

%--------------------------------------------------------------------------
    function AngleTree = RAY35_A(Pre,Int,SET)
        %Use ray-tracing to calculate 'angle tree' for textured interface.
        
        %INPUT
        %N1,N2      refractive index medium above/below interface
        %Int        struct with current interface info (coat.med, coat.thi, )
        %w          wavelength
        
        %NOTE: this A-part only uses the angle of refraction from the
        %precalculation table. Reflectance and transmittance data is used only in
        %part-B.
        

        %TODO: for clarity the ray-plotting was removed, but could be restored.
        
        nai=length(SET.bound_ang_interval)-1;                       %nr. of angular intervals
        depth = 5;                                 %ray-tracing depth (nr of bounces)
        d1 = 2^depth-1;                            %dim 1 of AngleTree determined by depth
        %AngleTree = zeros(d1,SET.ray_nr,2*nai); %initialize with zeros

        AngleTree = 88*ones(d1,SET.ray_nr,2*nai); %initialize with single bounce from bottom scenario
        AngleTree(1,:,:) = -0.5;                   %first bounce cos(60 degree) = 0.5, minus from bottom
        AngleTree(2,:,:) = 10.5;                   %reflected ray terminated at bottom
        AngleTree(3,:,:) = -9.5;                   %transmitted ray terminated at top
        
        for c1=1:nai            %for every angular interval
            %---rays incident from above---
            ray0 = OD0([SET.bound_ang_interval(c1),SET.bound_ang_interval(c1+1)],SET.ray_nr,Int.geo_p,-1,Int.xy); %init. rays
            AngleTree(:,:,c1) = TRACE(Pre,Int.geo_p,ray0,d1);                         %ray-tracing
        end
        if isfield(Int,'skip_rear') && Int.skip_rear == 0     %do rear side illumination only if ligh reaches there
            for c1=1:nai
                %---rays incident from below---
                ray0 = OD0([SET.bound_ang_interval(c1),SET.bound_ang_interval(c1+1)],SET.ray_nr,Int.geo_p,1,Int.xy); %init. rays
                AngleTree(:,:,nai+c1) = TRACE(Pre,Int.geo_p,ray0,d1);                    %ray-tracing
            end
        end
        %Angle tree data is ordered aoi 0-->90 (above), 0-->90 (below)
        %--------------------------------------------------------------------------
        function [ray]=OD0(theta,nr,G,dwn,xy)
            %set origin O0 and direction D0 of initial rays
            
            %INPUT
            %theta          inc. angle [rad]
            %nr             nr of rays
            %z              z-coordinate (floor or ceiling)
            %dwn            ray direction (-1=down, +1=up)
            %OUTPUT
            %O0, D0         origin, direction of initial ray (x,y,z)
            
            %inclination angle theta (measured from z-axis)
            Theta=theta(1)+(theta(2)-theta(1))*rand(nr,1); %random in interval
            Phi=2*pi*rand(nr,1);  %azimuth angle phi, random in [0,90]
            if dwn==1, z=G{1}.z(1); else, z=G{1}.z(2); end
            
            ray(:,1)=xy(1)*rand(nr,1);          %origin x
            ray(:,2)=xy(2)*rand(nr,1);          %origin y
            ray(:,3)=z*ones(nr,1);              %origin z
            ray(:,4)=sin(Theta).*cos(Phi);      %direction x
            ray(:,5)=sin(Theta).*sin(Phi);      %direction y
            ray(:,6)=dwn*cos(Theta);            %direction z
            ray(:,7)=1;                         %index of ray (parent = 1)
            
            %---for each ray determine its subdomain coordinate [bx,by]---
            for bx=1:size(G,1)
                for by=1:size(G,2)      %for every subdomain (in 2D array)
                    x=G{bx,by}.x;       %x-range of subdomain
                    y=G{bx,by}.y;       %y-range of subdomain
                    ix=ray(:,1)>=x(1) & ray(:,1)<=x(2)...  %rays in range
                        & ray(:,2)>=y(1) & ray(:,2)<=y(2);
                    
                    ray(ix,8)=bx;       %x-coord of subdomain containing ray
                    ray(ix,9)=by;      %y-coord of subdomain containing ray
                end
            end
            
        end
        
        %--------------------------------------------------------------------------
        function [AngleTree] = TRACE(Pre,G,ray0,d1)
            %tracing incident rays and recording angle of incidence upon every
            %bounce
            
            %INPUT
            %Pre        precalculation table (r,a,t vs angle of incidence)
            %G          geometry struct (pyramid points and surf. normals)
            %O0, D0     initial origin and direction of incident rays
            %plt        plot rays while tracing
            %OUTPUT
            %Fdi,Cdi    direction and intensity of rays reaching floor/ceiling
            %Ap, As     absorption in coating (p&s pol)
            %Lp, Ls     intensity of terminated rays
            %Nrx        nr. of intersections for each incident ray (vector)
            
            nrr=size(ray0,1);        %nr of rays
            AngleTree = 88*ones(d1,nrr);
            %Note:contains the cosine of angles, which is always between -1 and +1
            %These angles are used to lookup R,A,T. Any value outside -1 to +1
            %range will give NaN for R, A and T.
            %A value of 88 indicates sub-ray does not exist (terminated parent).
            %A value between -11 and -9 indicates a sub-ray terminated at ceiling.
            %A value between +9 and +11 indicates a sub-ray terminated at the floor.
            
            for ry=1:nrr             %for every incident ray
                rq=[ray0(ry,:);zeros(30,9)];      %ray queue (contains 1st ray + 30 empty slots)
                p=1;                 %pointer to ray queue
                %incident ray has index 1. Child_R has index 2*ix, Child_T has index 2*ix+1
                
                %reflected sub-ray is traced till death (transm sub-rays go in queue)
                while p>0                %while there are transmitted sub-rays in queue
                    
                    O =rq(p,1:3);        %origin
                    D =rq(p,4:6);        %direction
                    ix = rq(p,7);
                    b =rq(p,8:9);       %sub-domain box
                    p=p-1;               %reduce pointer
                    C=1;                 %C=1 ray is alive (C=0 ray is dead)
                    %disp('new ray')
                    
                    while C             %while the reflected sub-ray is still alive
                        %calculate next intersection point
                        [X,N,id]=intersexy(O,D,G{b(1),b(2)});
                        
                        switch id
                            case 3              %SUB-DOMAIN SIDE WALL
                                O=X;
                                if b(2)==1      %if also main-domain side wall
                                    D(2)=-D(2); %flip y-coordinate (reflect ray)
                                else
                                    b(2)=b(2)-1; %move to next sub-domain
                                end
                                
                            case 4              %SUB-DOMAIN SIDE WALL
                                O=X;
                                if b(1)==size(G,1)  %if also main-domain side wall
                                    D(1)=-D(1);     %flip x-coordinate (reflect ray)
                                else
                                    b(1)=b(1)+1;    %move to next sub-domain
                                end
                                
                            case 5              %SUB-DOMAIN SIDE WALL
                                O=X;
                                if b(2)==size(G,2)  %if also main-domain side wall
                                    D(2)=-D(2);     %flip y-coordinate (reflect ray)
                                else
                                    b(2)=b(2)+1;    %move to next sub-domain
                                end
                                
                            case 6              %SUB-DOMAIN SIDE WALL
                                O=X;
                                if b(1)==1      %if also main-domain side wall
                                    D(1)=-D(1); %flip x-coordinate (reflect ray)
                                else
                                    b(1)=b(1)-1;%move to next sub-domain
                                end
                                
                            case 2              %BOTTOM DOMAIN BOUNDARY
                                cai=-D*N';      %cos(angle of inc.)
                                AngleTree(ix,ry) = cai+10;          %store angle of incidence (cosine)
                                C=0;            %kill sub-ray
                                
                            case 1             %TOP DOMAIN BOUNDARY
                                cai=-D*N';      %cos(angle of inc.)
                                AngleTree(ix,ry) = cai-10;          %store angle of incidence (cosine)
                                C=0;           %kill sub-ray
                                
                            otherwise  %this is a 'real' intersection!
                                
                                O=X;        %intersection is new origin
                                
                                [Dr,Dt,cai,TIR]=reflexy(N,D,Pre); %direction of reflected and transmitted rays
                                AngleTree(ix,ry) = cai;       %store angle of incidence (cosine)
                                
                                if 2*ix < d1                  %if depth has not been exceeded
                                    if ~TIR                    %if no TIR
                                        p=p+1;                    %to next slot in queue...
                                        rq(p,:)=[X,Dt,2*ix+1,b];  %...transmitted ray is added
                                    end
                                    ix = 2*ix;                %index of reflected ray is double of parent
                                    D=Dr;                     %reflected ray continues
                                else
                                    C = 0;                    %kill sub-ray
                                    %disp('ray-terminated')
                                end
                        end
                    end
                end
            end
        end
        %..........................................................................
        function [x,n,id]=intersexy(O,D,G)
            %calculate intersection point (and plane id) for 1 ray
            
            %INPUT
            %O          origin of incident ray (x,y,z)
            %D          direction of incident ray (x,y,z)
            %G          geometry data (points and normals)
            %OUTPUT
            %X          intersection point (x,y,z)
            %N          local surface normal (x,y,z)
            %id         identifier of intersected plane
            
            
            OO=G.onr*O;                        %vectorize origin O
            
            %---clever triangle test---
            T=sum(G.N.*(G.P-OO),2)./(G.N*D');  %intersection times
            ct=T>G.small;                      %check for positive times
            %---x-check geometry points (7 to end)---
            X=OO+T*D;                          %intersection points
            dX=(X-G.P).*G.s;                   %xy-coord relative to corner1
            dS=dX(:,1)/G.mx+dX(:,2)/G.my;      %sum coord x+y
            cx=dX(:,1)>0 & dX(:,2)>0 & dS<1;   %check points
            %---box points(1 to 6) are x-checked separately---
            cx(1)=X(1,1)>G.x(1)&& X(1,1)<G.x(2)&& X(1,2)>G.y(1)&& X(1,2)<G.y(2);
            cx(2)=X(2,1)>G.x(1)&& X(2,1)<G.x(2)&& X(2,2)>G.y(1)&& X(2,2)<G.y(2);
            cx(3)=X(3,1)>G.x(1)&& X(3,1)<G.x(2)&& X(3,3)>G.z(1)&& X(3,3)<G.z(2);
            cx(4)=X(4,2)>G.y(1)&& X(4,2)<G.y(2)&& X(4,3)>G.z(1)&& X(4,3)<G.z(2);
            cx(5)=X(5,1)>G.x(1)&& X(5,1)<G.x(2)&& X(5,3)>G.z(1)&& X(5,3)<G.z(2);
            cx(6)=X(6,2)>G.y(1)&& X(6,2)<G.y(2)&& X(6,3)>G.z(1)&& X(6,3)<G.z(2);
            
            T(~ct|~cx)=Inf;     %assign inf time for not pass check t or x
            [~,id]=min(T);      %find closest intersection
            x=X(id,:);          %corresponding intersection point
            n=G.N(id,:);        %corresponding surface normal
        end
        %..........................................................................
        function [Dr,Dt,cai,TIR]=reflexy(N,D,Pre)
            %calculate direction of reflected and transmitted ray and r, a, t
            
            %INPUT
            %N          intersected plane normal (x,y,z)
            %D          incident ray direction (x,y,z)
            %Pre        precalculation table
            %OUTPUT
            %Dr,Dt      direction of reflected/transmitted ray (x,y,z)
            %RATp,RATs  reflectance, abs. , transm (p&s-pol)
            
            cai=-D*N';          %cos(angle of inc.)
            row=interp1x(Pre(:,1),Pre(:,2:3),cai);    %interpolate from table
            %---reflectance, absorptance, transmittance---
            
            %---direction of reflected ray---
            Dn=N*cai;           %surface normal component of D
            Dp=D+Dn;            %surface parallel component of D
            Dr=Dp+Dn;           %direction reflected ray
            %---direction of transmitted ray---
            cr=row(1);sr=row(2);%refraction cosine and sine from table
            TIR = 0;            %total internal reflection?
            if Dp==0             %if Dp=[0,0] --> not normable
                Dt=-cr*N*sign(cai); %normal incidence
            else
                Dt=-cr*N*sign(cai)+sr*normal(Dp);
                if cr < 1e-9, TIR = 1; end    %detect TIR when car = 0
            end
            Dt=normal(Dt);Dr=normal(Dr);  %normalize
            %..........................................................................
            function v=normal(V)
                %normalize vector (so x^2+y^2+z^2=1)
                v=V/sqrt(sum(V.^2));
            end
            %.........................................................................
            function yi=interp1x(xjo,yjo,xi)
                %super-mega fast interpolation (modified interp1)
                
                r=find(xjo>=xi,1,'last');
                u = (xi-xjo(r))/(xjo(r+1)-xjo(r));
                yi=yjo(r,:)+(yjo(r+1,:)-yjo(r,:))*u;
            end
        end
    end
%--------------------------------------------------------------------------
    function [Rpa,Rsa,Rpb,Rsb,Apa,Asa,Apb,Asb,Tpa,Tsa,Tpb,Tsb,logbook]= RAY35_B(AT,Pre,SET)
        %Calculate scatter matrices from local angles of incidence of all sub-rays
        %---INPUT---
        %AT         angle tree (cosine of angle for all sub-rays, -10 or +10 to
        %           indicate ray termination at floor or ceiling). Data is stored
        %           in 3D array (dim1: 2^depth-1, dim2: nr. rays, dim3: aoi)
        %Pre        precalculation table (R,A,T versus angle of incidence)
        %S          settings structure (contains angular interval settings)
        %---OUTPUT---
        %Rpa,Rsa,Rpb,Rsb    reflectance scatter matrix for p&s polarization and light incident from above & below
        %Apa,Asa,Apb,Asb    absorptance vector for p&s polarization and light incident from above & below
        %Tpa,Tsa,Tpb,Tsb    transmittance scatter matrix for p&s polarization and light incident from above & below
        
        nrc=(size(Pre,2)-7)/2;      %number of SUBcoatings
        nrr = size(AT,2);           %number of rays used
        D = size(AT,1);             %depth of ray-tracing (intensity relevant up to D) 2^n-1
        hD = (D+1)/2-1;             %half depth (R,A,T only relevant up to hD) 2^(n-1)-1
        
        %---use precalculation table to find R,A,T corresponding to every sub-ray angle---
        %Note: data stored in 4D array. Same dimensions as AT, but with 4th
        %dimension to accomodate R,A,T and dim1 only up to half depth.


%         RAT_Tree = zeros(hD,size(AT,2),size(AT,3),4+2*nrc);
% 
%         for rat = 1:(4+2*nrc)       %for every R,A,T (p&s and can be multi-coating)
%            RAT_Tree(:,:,:,rat) = interp1(Pre(:,1),Pre(:,rat+3),AT(1:hD,:,:));
%             %TODO: this is a time consuming step
%         end


        RAT_Tree = interp1(Pre(:,1),Pre(:,4:end),AT(1:hD,:,:));
        %took it outside the for loop, but limited speed gain
        
        %interp1 will return NaN for anything outside -1 to +1 range, such as
        %sub-rays terminated at ceiling/floor (their R,A,T is not relevant)
        
        nterm = AT(hD+1:end,:,:) >=-1 & AT(hD+1:end,:,:) <=1;      %non-terminated rays
        
        %---use R,A,T to calculate intensity of every sub-ray (and coating abs.)---
        %Note: data stored in same format as AngleTree.
        %Terminated rays have R,A,T = NaN, but DO have intensity! Intensity is
        %always calculated from parent R,A,T. This is done separately for p&s
        %polarized rays.
        
        Ip = zeros(size(AT));       Is = Ip;       %initialize ray intensity
        Abs_p = zeros(nrc,2*SET.nr_ang_interval); Abs_s = Abs_p; %absortance is 'accumulated'
        
        %sub-ray intensities are calculated iteratively (parent --> child -->...)
        %as would be done in a conventional ray-tracing simulation
        Ip(1,:,:) = 1;         %I_parent = 1 (need to divide by nr_rays in the end)
        Is(1,:,:) = 1;
        
        for d = 1:hD           %intensity is calculated until full depth (2d+1)
            %---intensity of RTchild rays---
            %p-polarized ray
            Ip(2*d  ,:,:) = Ip(d,:,:) .* RAT_Tree(d,:,:,1);       %I_Rchild = I_parent * R
            Ip(2*d+1,:,:) = Ip(d,:,:) .* RAT_Tree(d,:,:,2+nrc);   %I_Tchild = I_parent * T
            %s-polarized ray
            Is(2*d  ,:,:) = Is(d,:,:) .* RAT_Tree(d,:,:,3+nrc);   %I_Rchild = I_parent * R
            Is(2*d+1,:,:) = Is(d,:,:) .* RAT_Tree(d,:,:,4+2*nrc); %I_Tchild = I_parent * T
            %---absorptance in coatings---
            for c = 1:nrc               %for every coating (can be multiple)
                %p-polarized ray
                %for every ray at every aoi at current sub-ray level (d): Abs = I_parent * A
                Abs_d = squeeze(Ip(d,:,:) .* RAT_Tree(d,:,:,1+c));  %abs intensity at level d
                Abs_d(isnan(Abs_d)) = 0;    %non-existent sub-rays have NaN, set to 0
                Abs_p(c,:) = Abs_p(c,:) + sum (Abs_d);   %sum over all rays and accumulate
                %s-polarized ray
                Abs_d = squeeze(Is(d,:,:) .* RAT_Tree(d,:,:,3+nrc+c));  %abs intensity at level d
                Abs_d(isnan(Abs_d)) = 0;    %non-existent sub-rays have NaN, set to 0
                Abs_s(c,:) = Abs_s(c,:) + sum (Abs_d);   %sum over all rays and accumulate
            end
        end
        
        %accumulated absorptance becomes the absorptance vector (output)
        Apa = Abs_p(:,1:SET.nr_ang_interval)/nrr;
        Asa = Abs_s(:,1:SET.nr_ang_interval)/nrr;
        Apb = Abs_p(:,(SET.nr_ang_interval+1):2*SET.nr_ang_interval)/nrr;
        Asb = Abs_s(:,(SET.nr_ang_interval+1):2*SET.nr_ang_interval)/nrr;
        
        Ibhp = Ip(hD+1:end,:,:);                     %intensity in bottom half (no room for children due to limited depth)
        Ibhp(isnan(Ibhp))=0;                          %non-existent rays get intensity 0
        Lp = squeeze(sum(sum(Ibhp.*nterm)))/nrr;      %loss vector contains intensity of sub-rays not reached ceiling or floor
        Lpa = Lp(1:SET.nr_ang_interval)';                           %loss vector vs aoi from above
        Lpb = Lp((SET.nr_ang_interval+1):2*SET.nr_ang_interval)';                  %loss vector vs aoi from below
        
        Ibhs = Is(hD+1:end,:,:);                     %intensity in bottom half (no room for children due to limited depth)
        Ibhs(isnan(Ibhs))=0;                          %non-existent rays get intensity 0
        Ls = squeeze(sum(sum(Ibhs.*nterm)))/nrr;      %loss vector contains intensity of sub-rays not reached ceiling or floor
        Lsa = Ls(1:SET.nr_ang_interval)';                           %loss vector vs aoi from above
        Lsb = Ls((SET.nr_ang_interval+1):2*SET.nr_ang_interval)';                  %loss vector vs aoi from below
        
        %---reconstruct AID---
        Rpa = zeros(SET.nr_ang_interval); Rsa = Rpa; Rpb = Rpa; Rsb = Rpa;
        Tpa = Rpa;          Tsa = Rpa; Tpb = Rpa; Tsb = Rpa;   %initialize scatter matrices
        
        %---incidence from ABOVE---
        for aoi = 1:SET.nr_ang_interval              %for every angle of incidence
            
            ATaoi = AT(:,:,aoi);       %consider 1 angle of incidence (1 slice)
            Iaoi_p = Ip(:,:,aoi);      %consider 1 angle of incidence (1 slice)
            Iaoi_s = Is(:,:,aoi);
            
            ttx = ATaoi >= -11 & ATaoi <= -9;   %find rays reaching domain ceiling
            btx = ATaoi >= 0 & ATaoi <=11;      %find rays reaching domain floor
            
            cai_R = ATaoi(ttx)+10;     %cosine angle of exit to ceiling
            cai_T = ATaoi(btx)-10;     %cosine angle of exit to floor
            
            Iaoi_Rp = Iaoi_p(ttx);     %intensities of rays hitting ceiling
            Iaoi_Tp = Iaoi_p(btx);     %intensities of rays hitting floor
            Iaoi_Rs = Iaoi_s(ttx);     %intensities of rays hitting ceiling
            Iaoi_Ts = Iaoi_s(btx);     %intensities of rays hitting floor
            %indexing like this makes the format linear (a vector), and that's fine.
            
            for aoe = 1:SET.nr_ang_interval       %for every angle of exit interval
                %populate scatter matrices one element at the time
                ixr = cai_R < cos(SET.bound_ang_interval(aoe)) & cai_R > cos(SET.bound_ang_interval(aoe+1));   %find rays reflected in that interval
                Rpa(aoe,aoi) = sum(Iaoi_Rp(ixr))/nrr;                %sum their intensity
                Rsa(aoe,aoi) = sum(Iaoi_Rs(ixr))/nrr;                %sum their intensity
                
                ixt = cai_T < cos(SET.bound_ang_interval(aoe)) & cai_T > cos(SET.bound_ang_interval(aoe+1));   %find rays transmitted in that interval
                Tpa(aoe,aoi) = sum(Iaoi_Tp(ixt))/nrr;                %sum their intensity
                Tsa(aoe,aoi) = sum(Iaoi_Ts(ixt))/nrr;                %sum their intensity
            end
        end
        
        %---incidence from BELOW---
        for aoi = 1:SET.nr_ang_interval              %for every angle of incidence
            
            ATaoi = AT(:,:,SET.nr_ang_interval+aoi);       %consider 1 angle of incidence (1 slice)
            Iaoi_p = Ip(:,:,SET.nr_ang_interval+aoi);      %consider 1 angle of incidence (1 slice)
            Iaoi_s = Is(:,:,SET.nr_ang_interval+aoi);
            
            ttx = ATaoi >= -11 & ATaoi <= -9;   %find rays reaching domain ceiling
            btx = ATaoi >= 0 & ATaoi <=11;      %find rays reaching domain floor
            
            cai_T = ATaoi(ttx)+10;     %cosine angle of exit to ceiling
            cai_R = ATaoi(btx)-10;     %cosine angle of exit to floor
            
            %incidence from below, so reflected rays his floor, transmitted rays
            %hit ceiling!!!
            Iaoi_Rp = Iaoi_p(btx);     %intensities of rays hitting ceiling
            Iaoi_Tp = Iaoi_p(ttx);     %intensities of rays hitting floor
            Iaoi_Rs = Iaoi_s(btx);     %intensities of rays hitting ceiling
            Iaoi_Ts = Iaoi_s(ttx);     %intensities of rays hitting floor
            %indexing like this makes the format linear (a vector), and that's fine.
            
            for aoe = 1:SET.nr_ang_interval       %for every angle of exit interval
                %populate scatter matrices one element at the time
                ixr = cai_R < cos(SET.bound_ang_interval(aoe)) & cai_R > cos(SET.bound_ang_interval(aoe+1));   %find rays reflected in that interval
                Rpb(aoe,aoi) = sum(Iaoi_Rp(ixr))/nrr;                %sum their intensity
                Rsb(aoe,aoi) = sum(Iaoi_Rs(ixr))/nrr;                %sum their intensity
                
                ixt = cai_T < cos(SET.bound_ang_interval(aoe)) & cai_T > cos(SET.bound_ang_interval(aoe+1));   %find rays transmitted in that interval
                Tpb(aoe,aoi) = sum(Iaoi_Tp(ixt))/nrr;                %sum their intensity
                Tsb(aoe,aoi) = sum(Iaoi_Ts(ixt))/nrr;                %sum their intensity
            end
        end
        
        %---renormalization---
        %to correct for known loss L (terminated rays)
        %(not doing 'blind' renormalization to correct for unknowns)
        logbook = [];
        [Rpa,Apa,Tpa,logbook] = reno(Rpa,Apa,Tpa,Lpa,logbook);
        [Rsa,Asa,Tsa,logbook] = reno(Rsa,Asa,Tsa,Lsa,logbook);
        [Rpb,Apb,Tpb,logbook] = reno(Rpb,Apb,Tpb,Lpb,logbook);
        [Rsb,Asb,Tsb,logbook] = reno(Rsb,Asb,Tsb,Lsb,logbook);
        %Note: losses reduce with increasing depth. Also increasing the number of
        %rays helps (not sure why)
        
        %TO DO: check energy conservation in matrix and final R,A,T
        %TO DO: investigate how number of rays affects L
        %--------------------------------------------------------------------------
        %     function yi=interp1x(xjo,yjo,xi)
        %         %super-mega fast interpolation (modified interp1)
        %
        %         r=find(xjo>=xi,1,'last');
        %         u = (xi-xjo(r))/(xjo(r+1)-xjo(r));
        %         yi=yjo(r,:)+(yjo(r+1,:)-yjo(r,:))*u;
        %     end
        %--------------------------------------------------------------------------
        function [R,A,T,logbook]=reno(R,A,T,L,logbook)
            %renormalization (only very very very slight <1e-3)
           

            if max(L)>SET.ray_energy_conserv/SET.ray_nr
                %TODO: test division by ray_nr is ok
                logbook = [logbook,'Significant renormalization!'];
                %disp('Significant renormalization!');
                %TODO: also give energy loss value
            end
            R=R.*(ones(size(R,1),1)*(1./(1-L)));     %re-scaling R
            A=A.*(ones(size(A,1),1)*(1./(1-L)));     %re-scaling A
            T=T.*(ones(size(T,1),1)*(1./(1-L)));     %re-scaling T
        end
    end
%---------------------------------------------------------------------------
end