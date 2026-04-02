function [V,F,A,C,zenith,azimuth] = icohemisphere(rl,cut_off_half)
%Divide 'unit'-hemisphere into triangles and get the coordinates of triangle
%corners and centers. Triangles have a similar (but not identical size) and
%a total area of exactly 2pi. Useful for discretization of skydome. 

%By Rudi Santbergen, PVMD, TU Delft (2017), changed by Youri Blom

%INPUT
%rl             recursion level 1,2,3,... (nr. of triangles n = 10*4^rl)
%cut_off_half   indicates whether the bottom half should be removed.
%OUTPUT
%V              xyz-coord. of triangle vertex points [v x 3]
%F              triangle facets (see help on 'patch' function) [n x 3]
%A              triangle area [n x 1]
%C              xyz-coord. of triangle centers [n x 3]
%zenith         zenith angle of triangle centers [n x 1]
%azimuth        azimuth angle of triangle centers [n x 1]

%Strategy: based on icosphere algorithm, starting with icosahedron (20 
%triangles) each triangle is devided into 4 smaller ones, recusively.
%The icosphere obtained after 1 recursion (80 triangles) has mirror
%symmetry and is cut into hemisphere (40 triangles). This hemisphere is
%then refined recursively, multiplying the number of triangle by 4 at each
%recursion.

%---create icosahedron vertices and faces (20 triangles)---
t = (1+sqrt(5))/2;
V = [-1 t 0; 1  t  0; -1 -t  0; 1 -t  0; 0 -1  t; 0  1  t;
    0 -1 -t; 0  1 -t; t  0 -1; t  0  1; -t  0 -1; -t  0  1];

F = [1 12  6; 1  6  2; 1  2  8; 1  8 11; 1 11 12; 2  6 10; 6 12  5;
    12 11  3; 11  8  7; 8  2  9; 4 10  5; 4  5  3; 4  3  7;
    4  7  9; 4  9 10; 5 10  6; 3  5 12; 7  3 11; 9  7  8; 10  9  2];

%---rotate icosahedron to get mirror plane on xy plane---
V = rotx(V,atand(1/t));

%---refine once (this is needed to obtain mirror symmetry)---
[V,F] = refine(V,F);        %refine icosphere --> 80 trianlges

%---cut off bottom half (sphere becomes hemisphere--> 40 triangles)---
if cut_off_half; [V,F] = northern2(V,F); end

%---recusively refine until recursion level rl given as input---
for r = 1:rl-1              %for every recursion
    [V,F] = refine(V,F);    %refine icosphere
end

[A,C,zenith,azimuth] = postprocess(V,F); %get triangle center coordinates

%sum of triangles' area is always slightly less than 2 pi because it 
%is not a perfect sphere but approximation by triangles
A = A*2*pi/sum(A);          %rescale to make area exactly 2 pi
%......................................................................
    function V = rotx(V,ang)
        %rotate xyz coordinates in [v x 3] matrix V by angle ang x-axis
        
        M = [cosd(ang),-sind(ang);sind(ang),cosd(ang)]; %rotation matrix
        
        for c = 1:size(V,1)         %for every vertex             
            V(c,2:3) = V(c,2:3)*M;  %rotate
        end
    end
    %......................................................................
    function [V,F] = northern2(V,F)
        %takes a structure defined by vertices V and faces F and removes
        %all vertices with z<0. Faces that have one or more vertices with 
        %z<0 are removed entirely. Used to turn icosphere into icohemisphere
        
        xv = find(V(:,3)<0);      %find all southern vertices (with z<0)
        lv0 = size(V,1);          %original nr of vertices
        V(xv,:) = [];             %remove southern vertices
                
        xf = false(size(F,1),1);  %initialize southern facet index
        
        %indentify all facets that have at least one southern vertex
        for f = 1:size(F,1) %for every facet
            %is 1st, 2nd or 3rd corner of triangle a southern vertex?
            if ~isempty(find(xv==F(f,1),1)) ||~isempty(find(xv==F(f,2),1)) ||~isempty(find(xv==F(f,3),1))
                xf(f)=true;     %if yes, flag true
            end
        end
        
        F(xf,:) = [];            %remove flagged facets
        
        %Not done yet! Removal of southern vertices requires shifting of 
        %vertex numbers in facet matrix F!!!
        shift = 0;               %shift starts at 0
        for c = 1:lv0            %for every vertex
            if isempty(find(xv==c,1))%if vertex is NOT southern, shift index up according to shift history
                F(F==c) = c-shift;
            else
                shift = shift+1; %if vertex IS southern, index will be overwritten
            end
        end
    end
    %......................................................................
    function [Vnew,Fnew] = refine(V,F)
        %refines geometery by dividing every triangle into 4 smaller
        %triangles (same as icosphere algorithm)
        
        nrf = size(F,1);              %nr of faces
        Fnew = zeros(4*nrf,3);        %initialize new faces matrix
        fx = 1;                       %faces matrix index
        
        nrv = size(V,1);              %nr of vertices
        Vnew = [V; zeros(nrv,3)];     %initialize Vnew (half is old V)
        vx = nrv+1;                   %vertices matrix index
        
        S = zeros(nrv);               %side matrix
        
        for f = 1:nrf                 %for every facet
            f1 = F(f,1);              %vertex index 1
            f2 = F(f,2);              %vertex index 2
            f3 = F(f,3);              %vertex index 3
            
            s12 = S(f1,f2);           %side12
            if s12 == 0               %if side not yet divided
                v12 = (V(f1,:)+V(f2,:))/2; %calculate midpoint
                Vnew (vx,:) = v12;    %midpoint is new vertex
                S(f1,f2) = vx;        %remember that it is divided
                S(f2,f1) = vx;        %12 and 21 are the same
                f12 = vx;             %vertex index for F
                vx = vx+1;            %increment index
            else                      %if side already divided
                %midpoint already exists (no need to create)
                f12 = s12;            %vertex index for F
            end
            
            s23 = S(f2,f3);           %side23
            if s23 == 0               %if side not yet divided
                v23 = (V(f2,:)+V(f3,:))/2; %calculate midpoint
                Vnew (vx,:) = v23;    %midpoint is new vertex
                S(f2,f3) = vx;        %remember that it is divided
                S(f3,f2) = vx;        %23 and 32 are the same
                f23 = vx;             %vertex index for F
                vx = vx+1;            %increment index
            else                      %if side already divided
                %midpoint already exists (no need to create)
                f23 = s23;            %vertex index for F
            end
            
            s31 = S(f3,f1);           %side31
            if s31 == 0               %if side not yet divided
                v31 = (V(f3,:)+V(f1,:))/2; %calculate midpoint
                Vnew (vx,:) = v31;    %midpoint is new vertex
                S(f3,f1) = vx;        %remember that it is divided
                S(f1,f3) = vx;        %31 and 13 are the same
                f31 = vx;             %vertex index for F
                vx = vx+1;            %increment index
            else                      %if side already divided
                %midpoint already exists (no need to create)
                f31 = s31;            %vertex index for F
            end
            
            %add 4 new facets
            Fnew(fx,:)   = [f1,f12,f31];
            Fnew(fx+1,:) = [f2,f23,f12];
            Fnew(fx+2,:) = [f3,f31,f23];
            Fnew(fx+3,:) = [f12,f23,f31];
            fx = fx+4;                %increment
            
            Vnew = normal(Vnew);        %renormal after each refinement
            %refining only after final refinement gives less uniform
            %surface area
        end
        %..................................................................
        function V = normal(V)
            %---normalize to unit radius---
            n = sqrt(V(:,1).^2+V(:,2).^2+V(:,3).^2);
            V = V./(n*[1,1,1]);
        end
    end
    %......................................................................
    function [A,C,zenith,azimuth] = postprocess(V,F)
        %Calculate face center xyz and face area
        
        V0 = V(F(:,1),:);          %first vertex
        V1 = V(F(:,2),:);          %first vertex' one side neighbor
        V2 = V(F(:,3),:);          %first vertex' other side neighbor
        
        E1 = V1-V0;                %edge vector1
        E2 = V2-V0;                %edge vector2
        
        %---edge vector cross-product gives surface normal vector---
        Nx = E1(:,2).*E2(:,3)-E1(:,3).*E2(:,2);     %x-component
        Ny = E1(:,3).*E2(:,1)-E1(:,1).*E2(:,3);     %y-component
        Nz = E1(:,1).*E2(:,2)-E1(:,2).*E2(:,1);     %z-component
        %sign convention such that normal points to RIGHT HAND RULE side
        
        A = sqrt(Nx.^2+Ny.^2+Nz.^2)/2;     %area of triangle
        
        C = (V0+V1+V2)/3;                  %vertex center xyz-coord
        
        zenith = acosd(C(:,3));            %zenith angle 0,90 [degree]
        azimuth = atan2d(-C(:,1),-C(:,2));  %azimuth angle -180,180 [degree]
        %swapped x-->-y and y-->-x such that 
        %0   deg is south
        %90  deg is west
        %180 deg is north
        %270 deg is east
        
    end
end
%--------------------------------------------------------------------------