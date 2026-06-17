function [W,G]= pRAY02(Z,xy)
%Geometry data (points, normals, etc.) for ray-tracing. For speed
%the domain is devided into subdomains and the data is stored in a
%2D-struct.

%INPUT
%Z          height data matrix [um]
%xy         size heigth data x&y [um]
%OUTPUT
%G          2D-struct containing data for each subdomain separately
%W          1x1 struct containing data of whole domain

%if numel(xy)==1, xy = xy*[1,1];end  %if single xy, assume square x = y.

zbs=1e-9;                       %box z-spacing
zf=min(min(Z))-zbs;             %z coord floor (same for every box)
zc=max(max(Z))+zbs;             %z coord ceiling (same for every box)

[npx,npy]=size(Z);              %nr of points x&y direction
nsx=floor(sqrt(npx));           %nr of subdomains x (empirical)
nsy=floor(sqrt(npy));           %nr of subdomains y (empirical)
G(nsx,nsy)={struct([])};        %initialization
%TO DO: check for no-mans land in between the domains?
%subdomains are numbered by x=1,2,... and y=1,2,3...
for by=1:nsy                    %for every y-coordinate
    iy1=1+round((by-1)*(npy-1)/nsy);    %subdomain first point y
    iy2=1+round( by   *(npy-1)/nsy);    %subdomain end point y
    y1=(iy1-1)*xy(2)/(npy-1);           %subdomain start y [m]
    y2=(iy2-1)*xy(2)/(npy-1);           %subdomain end y [m]
    
    for bx=1:nsx                %for every x-coordinate
        ix1=1+round((bx-1)*(npx-1)/nsx);  %subdomain first point x
        ix2=1+round( bx   *(npx-1)/nsx);  %subdomain end point y
        x1=(ix1-1)*xy(1)/(npx-1);         %subdomain start x [m]
        x2=(ix2-1)*xy(1)/(npx-1);         %subdomain end x [m]
        %each time fill 2D-struct G with subdomain data
        G{bx,by}=METRY(Z(ix1:ix2,iy1:iy2),[x1,x2],[y1,y2],[zf,zc]);
    end
end
%once fill struct W with whole domain data (for plotting)
W=METRY(Z,[0,xy(1)],[0,xy(2)],[zf,zc]);

% if SET.plot_fig>1
%     %---plot geometry---
%     n = 70+i;               %figure nr.
%     figure(n)
%     clf(n)            %clear figure
%     %trisurf(W.tri,W.X,W.Y,W.Z,'FaceColor','red','FaceAlpha',0.3) %textured surface
%     h = trisurf(W.tri,W.X,W.Y,W.Z,'FaceLighting','gouraud'); %textured surface
%     shading flat
%     h.FaceColor = 'red';        %shading resets color, so done after
%     camlight('right')
%     patch('Vertices',W.Vbx,'Faces',W.Fbx,'FaceVertexCData',...
%         [0,1,0],'FaceColor','flat','FaceAlpha',0.1)              %transparent box
%     axis equal; xlim([W.x(1),W.x(2)]); ylim([W.y(1),W.y(2)]);
%     zlim([W.z(1),W.z(2)]);
%     %set(n,'Position',[200 200 500 400])
%     view(30,20)
%     title(['Height map interface ',num2str(i)])
%     xlabel('x [\mum]')
%     ylabel('y [\mum]')
%     zlabel('z [\mum]')
%     drawnow
% end
%TO DO: make plot also for WAVE model? No! This is to check triangulation.
%..............................................................
    function [g]=METRY(Z,x,y,z)
        %prepare geometry data as required for ray-tracing
        %INPUT
        %Z          height data matrix [um]
        %x,y,z      bounding box coordinates x=[x1,x2] [um]
        
        %---bounding box points and normals---
        cor1=[x(1),y(1),z(1)];cor2=[x(2),y(2),z(2)];   %opposite corners
        Pb=[cor2;cor1;cor1;cor2;cor2;cor1];            %6 box points
        Nb=[0,0,-1;0,0,1;0,1,0;-1,0,0;0,-1,0;1, 0, 0]; %6 box normals
        
        %---height data points and normals---
        [nrx,nry]=size(Z);            %nr of elements x&y
        g.mx=(x(2)-x(1))/(nrx-1);       %mesh size x
        g.my=(y(2)-y(1))/(nry-1);       %mesh size y
        [X,Y]=meshgrid(x(1):g.mx:x(2),y(1):g.my:y(2)); %create X,Y grid
        X=X.'; Y=Y.';                 %transpose
        
        %---make lx3 matrix containing every triangle, each row
        %containing linear index of 1st, 2nd and 3rd corner
        %---south west triangles---
        tri1=zeros(1,(nrx-1)*(nry-1));  %initialize
        for c=1:nry-1                   %for every triangle
            tri1((c-1)*(nrx-1)+(1:nrx-1))=(c-1)*nrx+(1:nrx-1); %1st
        end
        tri2=tri1+1;                    %2nd corner linear index
        tri3=tri1+nrx;                  %3rd corner linear index
        
        %---north east triangles (will have opposite normals!)
        tri6=zeros(1,(nrx-1)*(nry-1));  %initialize
        for c=1:nry-1                   %for every triangle
            tri6((c-1)*(nrx-1)+(1:nrx-1))=(c-1)*nrx+(2:nrx); %1st
        end
        tri5=tri6+nrx-1;                %2nd corner linear index
        tri4=tri6+nrx;                  %3rd corner linear index
        
        g.tri=[tri1',tri2',tri3';tri4',tri5',tri6']; %lin. index matrix
        %---convert linear index to list of coord--
        Xt=X(g.tri);                      %lx3 list of x-coord
        Yt=Y(g.tri);                      %lx3 list of y-coord
        Zt=Z(g.tri);                      %lx3 list of z-coord
        %---
        P=[Xt(:,1),Yt(:,1),Zt(:,1)];    %triangle POINTS: all 1st corners
        g.P=[Pb;P];                     %add box POINTS
        %---v=corner 1 minus corner 2---
        vx=Xt(:,1)-Xt(:,2);             %x-coord
        vy=Yt(:,1)-Yt(:,2);             %y-coord
        vz=Zt(:,1)-Zt(:,2);             %z-coord
        %---u=corner 1 minus corner 3---
        ux=Xt(:,1)-Xt(:,3);             %x-coord
        uy=Yt(:,1)-Yt(:,3);             %y-coord
        uz=Zt(:,1)-Zt(:,3);             %z-coord
        %---triangle NORMAL N = v x u (cross product)---
        N(:,1)=vy.*uz-vz.*uy;           %x-coord
        N(:,2)=vz.*ux-vx.*uz;           %y-coord
        N(:,3)=vx.*uy-vy.*ux;           %z-coord
        l=sqrt(sum(N.^2,2))*[1,1,1]; N=N./l; %normalize
        g.N=[Nb;N];       %add box normals
        
        %===additional parameters===
        %---for plotting---
        g.X=X;g.Y=Y;g.Z=Z;      %texture coordinates X,Y,Z
        g.Vbx=[x(1),y(1),z(1);x(2),y(1),z(1);x(2),y(2),z(1);...
            x(1),y(2),z(1);x(1),y(1),z(2);x(2),y(1),z(2);...
            x(2),y(2),z(2);x(1),y(2),z(2)];     %box vertices
        g.Fbx=[1,2,6,5;2,3,7,6;3,4,8,7;4,1,5,8;...
            1,2,3,4;5,6,7,8];                   %box faces
        %---
        g.x=x;g.y=y;g.z=z;      %input x,y,z bounding box [start,end]
        g.small=1e-12;          %small nr
        nr=size(g.N,1);         %nr. of triangles
        g.onr=ones(nr,1);       %row of ones
        %---for intersection test (distinguish sw and ne triangle)---
        ne=nr/2+4:nr;           %north-east triangle index
        g.s=ones(nr,3);         %3 rows of ones
        g.s(ne,:)=-1;           %make north-east part of 3rd row -1
        %---
    end
end