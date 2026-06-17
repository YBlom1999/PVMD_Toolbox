function [Rpa,Rsa,Rpb,Rsb,Apa,Asa,Apb,Asb,Tpa,Tsa,Tpb,Tsb]=FLAT17(N1,N2,coat,SET)
%Calculate R,A,T matrices for flat interface. Both for light coming above
%and below and both for p and s polarization.
%INPUT:
%N1,N2      refractive indices on either side of the flat interface
%coat       coating struct for interface
%S          settings struct (mostly from GP4_settings.m)
%OUTPUT:
%Rpa,Rsa,Rpb,Rsb,Apa,Asa,Apb,Asb,Tpa,Tsa,Tpb,Tsb    R,A,T matrices for flat interface

nrc=length(coat);                 %number of coatings

%---prepare refractive index and thickness vectors---
NN=zeros(2+nrc,1);                %refractive index vector (preallocate)
%tt=NN;                            %thickness vector (preallocate)
tt = cell(2+nrc,1);                %thickness cell array (initialize)

NN(1)=N1;                         %incident medium N
%tt(1)=inf;                        %incident medium t

for c=1:nrc                       %for every coating
    NN(1+c) = coat(c).N(SET.w);       %coating N
    %tt(1+c)=coat(c).thi(end);     %coating t
    tt{1+c} = coat(c).thi;     %coating t
end

NN(2+nrc)=N2;                     %outgoing medium N
%tt(2+nrc)=inf;                    %outgoing medium t

%---calculate R,A,T matrices---
[Rsa,Rpa,Asa,Apa,Tsa,Tpa]=augustin(NN,tt,SET);                  %light coming in from top
[Rsb,Rpb,Asb,Apb,Tsb,Tpb]=augustin(flipud(NN),flipud(tt),SET);  %light coming in from bottom (fliped NN and tt vectors)

Apb=flipud(Apb);Asb=flipud(Asb);    %coating is flipped!!!
%==========================================================================

    function [Rs,Rp,As,Ap,Ts,Tp]=augustin(NN,tt,SET)
        %Calculates R,A,T matrices for flat interface that may contain a
        %multilayer coating. Fresnel02.m calculates r,a,t values that go 
        %into the matrices.
        %INPUT
        %NN         refractive index vector
        %tt         thickness vector
        %aib        angular interval boundaries (assumed uniform!)
        %wav        wavelength
        %OUTPUT
        %Rs,Rp      reflection matrix (p&s pol)     [nai x nai] 
        %Ap,As      absorpion matrix (p&s pol)      [nrc x nai]
        %Tp,Ts      transmittance matrix (p&s pol)  [nai x nai] 
        
        nrs=10;                         %nr of sub-rays per interval
        nai=length(SET.bound_ang_interval)-1;              %number of angular intervals
        ad=90/(nai*nrs);                %angular distance between rays (deg)
        rai=((ad/2):ad:(90-ad/2))*pi/180;   %sub-ray angles (deg)
        
        %calculate r,a,t and angle of refraction for every sub-ray (p&s)
        [RATp,RATs,arr]=fresnel06(NN,tt,rai.',SET.wavelength_um(SET.w));
        
        %---preallocate R, A and T matrices---
        Rs=zeros(nai);      %reflection matrix s-pol           
        Rp=Rs;              %reflection matrix p-pol
        Ts=Rs;              %transmission matrix s-pol 
        Tp=Rs;              %transmission matrix p-pol
        
        As=zeros(size(RATs,1)-2,nai); %absorption matrix s-pol
        Ap=As;                         %absorption matrix s-pol
        
        
        
        for c1=1:nai                               %for each incident interval...
            ix=rai>SET.bound_ang_interval(c1)&rai<SET.bound_ang_interval(c1+1);          %index of rays incident in interval c1
            six = sum(ix);                         %nr of sub-rays in incident interval (for uniform aib six = nrs)  

            As(:,c1)=mean(RATs(2:end-1,ix),2);     %A-matrix s-pol (one vector per coating layer)
            Ap(:,c1)=mean(RATp(2:end-1,ix),2);     %A-matrix p-pol (one vector per coating layer)
            
            Rs(c1,c1)=mean(RATs(1,ix),2);          %R-matrix s-pol (diagonal matrix)
            Rp(c1,c1)=mean(RATp(1,ix),2);          %R-matrix p-pol (diagonal matrix)
            
            ts1=RATs(end,ix);                      %transmittance s-pol
            tp1=RATp(end,ix);                      %transmittance p-pol
            rar1=arr(ix);                          %angles of refraction 
            str=find(SET.bound_ang_interval<=rar1(1),1,'last');       %start c2 interval
            stp=find(SET.bound_ang_interval>=rar1(end),1,'first')-1;  %stop c2 interval
            for c2=str:stp                         %for each RELEVANT outgoing interval...
                %calculate irradiance received from sub-rays incident in c1 going to c2
                Ts(c2,c1)=sum(ts1(rar1>SET.bound_ang_interval(c2)&rar1<SET.bound_ang_interval(c2+1)))/six; %T-matrix s-pol
                Tp(c2,c1)=sum(tp1(rar1>SET.bound_ang_interval(c2)&rar1<SET.bound_ang_interval(c2+1)))/six; %T-matrix p-pol
            end
        end
    end

end