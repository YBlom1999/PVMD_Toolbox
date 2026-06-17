function Ni=read_nk(MaterialName,wavelengths)
%Imports refractive index from file and interpolates on specified wavelength

%--INPUT--
%material        material e.g. 'c-Si' (should be file name or constant index)
%wavelengths     wavelength vector for interpolation [um]

%--OUTPUT--
%Ni              complex refractive index vector interpolated at wavelength_i

flp = 0;
nkFilePath = fullfile('materials',[MaterialName,'.xlsx']);   % xlsx materials folder
if exist(nkFilePath,'file')           %===if file exists

    T = readtable(nkFilePath,"Sheet",'nk');

    %---read wavelength (1st column) and convert to um---
    colunn1 = T{:,1};                           %get wavelength data
    wl_unit = T.Properties.VariableNames{1};    %get wavelength units (um/nm/m/eV)

    switch wl_unit
        case 'um'               %units in micron
            wl_um = colunn1;
        case 'nm'               %units in nanometre
            wl_um = colunn1/1e3;
        case 'm'                %units in metre
            wl_um = colunn1*1e6;
        case 'eV'               %units in electron-volt
            wl_um = flipud(1.2398./colunn1);
            flp = 1;            %remember to flip nk-data as well
        otherwise
            error(['Wavelength units in ',nkFilePath,'?']);
    end

    %---read nk (2nd and 3rd column) and convert alpha/epsilon---
    column2 = T{:,2};                           %get n data
    column3 = T{:,3};                           %get k data
    k_unit = T.Properties.VariableNames{3};     %get k units (k/alpha/epsilon)

    switch k_unit       %column 3 header
        case 'k'                 %in k format
            n = column2;         %real
            k = column3;         %imaginary
            N = n + 1i*k;        %complex refractive index

        case 'alpha'             %in alpha format, abs. coef, assuming units [cm-1]
            n = column2;         %real
            k = column3.*wl_um/(4e4*pi);  %imaginary (alpha to k, and cm to um)
            N = n + 1i*k;        %complex refractive index

        case 'e2'                %complex dielectric function
            epsilon1 = column2;        
            epsilon2 = column3;
            epsilon = epsilon1+1i*epsilon2;   %assume col2 = real, col3 = imag
            N = sqrt(epsilon);                %complex index
        otherwise
            error(['nk units in ',nkFilePath,'?']);
    end
    
    if flp, N = flipud(N); end   %flip nk data if eV

    if wavelengths(1)<wl_um(1) || wavelengths(end) > wl_um(end) 
        disp(['Wavelengths out of "',MaterialName,'" data range. Extrapolating!'])
    end
    
    %---wavelength INTERPOLATION---
    % interp1 has extrapolation same as interpolation method (linear) or 
    % predefined fixed value. Doing 'nearest' extrapolation manually. 

    less = wavelengths < wl_um(1);     %any wl less than minimum in file
    more = wavelengths > wl_um(end);   %any wl more than maximum in file
    fine = ~less & ~more;               %wl inside file data range

    Ni(less) = N(1);          %if any less, assign first N value
    Ni(fine) = interp1(wl_um,N,wavelengths(fine),"linear");    %linear interpolation
    Ni(more) = N(end);        %if any more, assign final N value
    %---
    
else            %===if not a file, it could be a (complex) number
    
    Nfix = str2double(MaterialName);        %convert string to number

    if isnan(Nfix)                      %if not a number
        error([MaterialName,' is not a file or a complex number.']);
    else
        %if (complex) number, assign it to refractive index (constant for
        %all wavelengths)
        Ni = Nfix*ones(size(wavelengths));    
    end
end

end