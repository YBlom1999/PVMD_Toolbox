%EXAMPLE4: varying layer thickness
%run example 3 first so that scatter matrices are stored in 'Int.SM'

%Do not clear 'Int' to preserve scatter matrices
S.plot = 1;                     %don't plot height map, only abs.
T = 50:25:500;                  %specify wafer thicknesses [um]
J = zeros(size(T));             %initialize

for t = 1:length(T)             %for every thickness
    disp(['thickness: ',num2str(T(t)),' um']) %display
    Lay(2).thi = T(t);          %adjust layer thickness
    [Lay,Int] = GENPRO4(Lay,Int,S); %re-run simulation
    J(t) = Lay(2).cur;           %store current in vector
    drawnow
end

%plot current versus thickness
figure(1)                        
plot(T,J,'or')                       
xlabel('c-Si wafer thickness [\mum]')
ylabel('implied photocurrent [mA/cm^2]')