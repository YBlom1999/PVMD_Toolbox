function [RATp,RATs,roar]=fresnel06(N,thck,phi0,w)
%Calculates reflectance, absorptance and transmittance as a function of
%incident angle for FLAT and COHERENT multilayer structure. Uses
%net-radiation method applied to complex electric field strength [1].
% ---INPUT---
% N         refractive index vector (at single wavelength w) [-]
% t         layer thickness vector/cell array [um]
% phi0      angle of incidence vector (RADIANS!)
% w         single wavelength [um]

% ---OUTPUT---
% RATp,RATs refl, abs., transm. table for p&s polarized light [-]
%           every column corresponds to an incident angle, every row to a
%           (sub)layer (except 1st and final columns, which contain R and T)
% roar      angle of refraction (real) [rad]


[tv,tc] = thickness(thck);    %thickness cell and vector

%---calculation of real angle of refraction---
s=(N(1)/N(end))*sin(phi0);    %sin(angle of refraction)
c=sqrt(1-s.^2);               %cos(angle of refraction)
roar=atan(sin(phi0)./real(N(end)*c/N(1))); %real angle of refraction
%note: roar calculation is independent from RAT calculation

%---individual interface (r,t) and layer (tau), all complex numbers---
[rt_p,rt_s,tau,Np,Ns] = optics(N,tv,phi0,w);

%---solve equations to get fluxes for every angle of inc---
Qp = zeros(4*(length(N)-1),length(phi0)); Qs = Qp; %initialize fluxes matrix

for a=1:length(phi0)                         %for every angle of incidence
    Qp(:,a)=solvit(rt_p(:,:,a),tau(:,a));    %solve set of eq. to obtain fluxes
    Qs(:,a)=solvit(rt_s(:,:,a),tau(:,a));    %for p&s polarization
end

%---get overall R,A,T from fluxes
RATp = powrat(Qp,Np,tc,Ns/w);    %<--Ns/w (this is NOT a mistake)
RATs = powrat(Qs,Ns,tc,Ns/w);    %<--Ns/w

%==========================================================================

    function [tv,tc] = thickness(thck)
        % two types of input can be given
        % t = [inf,0.090,0.040,inf]; VECTOR means 90 & 40 nm coating, no profiling
        % so gives 2 absorptance values (backwards compattible with fresnel02)
        % t = {inf,0:0.005:0.090,0:0.005:0.040,inf} CELL means 90 & 40 nm coating,
        % both with 5 nm profiling (giving 18 + 8 absorptance values)
        
        if iscell(thck)             %if thck is a cell, construct the vector tv
            tc = thck;              %cell is input
            tv = Inf(size(thck));   %initialize tv
            for l = 2:length(tv)-1
                tv(l) = tc{l}(end); %fill tv
                if numel(tc{l})==1, tc{l} = [0,tc{l}]; end
                %if the cell contains a single thickness value, add 0!!!
            end 
        else                        %if thck is a vector, construct the cell
            tv = thck;              %vector is input
            tc = cell(size(thck));  %initialize tc
            for l = 2:length(tc)-1, tc{l} = [0,tv(l)]; end %fill tc
        end
    end
%--------------------------------------------------------------------------
    function [rt_p,rt_s,tau,Np,Ns]=optics(N,t,phi0,w)
        %Calculate interface r, t and layer tau for all wavelengths and all
        %angles of incidence simultaneously
        % INPUT:
        % N         refractive index stack (at single wavelength) [-]
        % t         layer thickness stack [um] VECTOR
        % w         single wavelength [um]
        % phi0      angle of incidence [rad]
        % OUTPUT:
        % rt_p,rt_s rt-coef. (ra,rb,ta,tb for each interface and angle)[-]
        % tau       layer transmittance (for each interface and angle) [-]
        % Np,Ns     'tilted' N (see [2]) [-]
        % Note: all output are complex numbers
        
        N(1)=real(N(1));             %make incident medium real
        
        O1 = ones(length(N),1);      %dim1: nr. of layers
        O2 = ones(1,length(phi0));   %dim2: nr. of angles
        N2 = (N.^2)*O2;              %square and make N [dim1 x dim2]
        
        %---calculation of r and t---
        %N multiplied or divided by cosine of complex angle or refraction,
        %using Snell's law (but in fast, less familiar form), gives...
        Ns = sqrt(N2-(N(1)*O1*sin(phi0')).^2);   %'tilted' N for s-pol. [1]
        Np = N2./Ns;                           %'tilted' N for p-pol.
        
        rt_p = fresnel(Np);            %complex fresnel coef. ...
        rt_s = fresnel(Ns);            %for p&s polarisation
        
        %---calculation of tau---
        tau = exp(1i*2*pi*Ns.*(t*O2)/w);  %complex transmittance
        %the N used for optical thickness calculation N*t happens to be the
        %same as Ns (tilted for s-pol). Just coincidence.
        %..........................................................................
        function [rt]=fresnel(N)
            %fresnel equations compattible with [2]
            %3D-array: dim1:ra,ta,rb,tb, dim2:interface, dim3:angle
            
            n1=N(1:end-1,:);n2=N(2:end,:);
            ra=(n1-n2)./(n1+n2);
            rt(1,:,:)=ra; rt(2,:,:)=1+ra; %r&t above
            rt(3,:,:)=-ra;rt(4,:,:)=1-ra; %r&t below
        end
    end
%--------------------------------------------------------------------------
    function Q=solvit(rt,tau)
        %Coherent net-radiation method solves set of linear equations to
        %get net-radiation fluxes Q (= E-field strengths). Done for each
        %angle separately.
        
        I = size(rt,2);             %nr of interfaces
        M=diag(ones(4*I,1));        %start with diagonal matrix
        for i=1:I                   %for each interface (and layer below)
            %add r & t coefficients to matrix M
            x=((i-1)*4);
            M(x+2,x+1)=-rt(1,i);   %ra (reflectance above)
            M(x+2,x+3)=-rt(4,i);   %tb (transmittance below)
            M(x+4,x+1)=-rt(2,i);   %ta (transmittance above)
            M(x+4,x+3)=-rt(3,i);   %rb (reflectance below)
            if i~=I                %if not final interface
                M(x+3,x+6)=-tau(i+1);       %add tau (layer transm.) coef
                M(x+5,x+4)=-tau(i+1);
            end
        end
        C=zeros(4*I,1);C(1)=1;     %flux-vector (light incident from top)
        Q=M\C;              %solve set of eq. by matrix right devision
    end
%--------------------------------------------------------------------------
    function [RAT]=powrat(Q,N,tc,Nw)
        %extract reflectance, absorptance and transmittance from fluxes Q
        %this is the only part where the depth profile needs to be taken into
        %account. Done for all angles simultaneously.
        
        L = size(N,1);                    %nr of layers
        for l = L:-1:1, rows(l) = max(1,length(tc{l})-1); end %nr rows in table
        
        RAT = zeros(sum(rows),size(N,2)); %initialize RAT-table
        
        RAT(1,:) = abs(Q(2,:)).^2;        %reflectance in first column
        RAT(end,:) = real(N(L,:)./N(1,:)).*abs(Q(end,:)).^2; %transm. in final column
        
        o2 = ones(1,size(N,2));
        for l=2:L-1                       %for each layer (excl. first & final)
            o1 = ones(length(tc{l}),1);
            z1 = tc{l}'*o2;                %distances to coating front
            z2 = tc{l}(end)-z1;            %distances to coating rear
            x=((l-2)*4);
            E1 = exp(1i*2*pi*z1.*(o1*Nw(l,:))).*(o1*Q(x+4,:)); %E1-field along depth
            E2 = exp(1i*2*pi*z2.*(o1*Nw(l,:))).*(o1*Q(x+6,:)); %E2-field along depth
            E = E1 + E2;                  %total E-field
            H= (o1*(N(l,:)./N(1,:))).*(E1-E2);     %total H-field (down minus up)
            S=real(E.*conj(H));                    %pointing vector (power flux)
            RAT((sum(rows(1:l-1))+1):sum(rows(1:l)),:) = S(1:end-1,:)-S(2:end,:);
            %absorptance is power flux at interval top minus bottom
        end
        
        %check whether R, A and T are >0 and <1, if not give warning
        if min(min(RAT))<-1e-4 || max(max(RAT))>1+1e-4
            warning('MATLAB:wN','R, A or T is <0 or >1')
        end
        %check whether R+A+T=1, if not give warning.
        if sum(abs(sum(RAT)-1))>1e-3
            warning('MATLAB:wE','R + A + T ~=1');
        end
    end
%--------------------------------------------------------------------------
end

%[1]    R. Santbergen et al., Optics Express 21 (102) A262-267.
%[2]    H.A. MacLeod, Thin-Film Optical Filters, Third Edition (2001), ch2