function Ni=readnk(str,wli)
%Imports complex refractive index from file or database and interpolates
%value on specified wavelength
%INPUT
%str        material string (name or constant index)
%wli        wavelengths for interpolation [um]
%OUTPUT
%Ni         complex refractive index at wli wavelenghts

%THIS PART OF THE CODE IS 'OPEN-SOURCE' TO ALLOW USER EDITING

%TODO: importdata does not handle double tabs well

flp = 0;
file = ['nk/',str,'.nk'];       %nk file name
if exist(file,'file')           %IF FILE EXISTS
    %TODO: should be able to read standard ASA nk-files
    D = importdata(file);       %import data from file
    %---wavelength---
    switch D.colheaders{1}      %column 1 header
        case 'um'               %units in micron
            wl_um = D.data(:,1);
        case 'nm'               %units in nanometre
            wl_um = D.data(:,1)/1e3;
        case 'm'                %units in metre
            wl_um = D.data(:,1)*1e6;
        case 'eV'               %units in electron-volt
            wl_um = flipud(1.2398./D.data(:,1));
            flp = 1;            %remember to flip nk-data as well
        otherwise
            error(['Wavelength units in ',file,'?']);
    end
    %---nk---
    switch D.colheaders{3}       %column 3 header
        case 'k'                 %k-data
            n = D.data(:,2);
            k = D.data(:,3);
        case 'alpha'             %abs. coef [cm-1]!!!
            n = D.data(:,2);
            k = D.data(:,3).*wl_um/(4e4*pi);
        case 'e2'                %complex dielectric function
            eps = D.data(:,2)+1i*D.data(:,3);
            N = sqrt(eps);
            n = real(N);
            k = imag(N);
        otherwise
            error(['nk units in ',file,'?']);
    end
    
    if flp, n = flipud(n); k = flipud(k); end   %flip nk data if eV
    
    Ni = interpol_fun(wl_um,n,wl_um,k,wli,str); %interpolate to wli
    
else            %IF NO FILE, IS IT A (COMPLEX) NUMBER?
    
    Nfix = str2double(str);
    if isnan(Nfix)
        error([str,' is not a file or a complex number.']);
    else
        Ni = Nfix*ones(size(wli));       %complex refractive index
    end
end
%----------------------------------------------------------------------
    function Ni=interpol_fun(wln_um,n,wlk_um,k,wli,str)
        %actual interpolation function
        %INPUT
        %wln_um, wlk_um     wavelength vector for nk data
        %n, k               corresponding nk data
        
        ni=interp1(wln_um,n,wli,'linear','extrap');     %real part
        ki=interp1(wlk_um,k,wli,'linear','extrap');     %imaginary part
        Ni=ni+1i*ki;                            %complex refractive indes
        %Hint: TB could cause error below     
        if wli(1)<wln_um(1) || wli(end) > wln_um(end) || ...
                wli(1) < wlk_um(1) || wli(end) > wlk_um(end)
            disp(['Wavelengths out of "',str,'" data range. Extrapolating!'])
        end
    end
end