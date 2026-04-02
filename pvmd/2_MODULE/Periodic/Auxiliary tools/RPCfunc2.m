function [NrRays] = RPCfunc2(sd,l0,w0,L0,W0)
%---INPUT---
%sd     Desired standard deviation
%l0,w0  Desired cell size
%L0,W0  Desired unit cell size
%---OUTPUT---
%NrRays Suggested number of rays

%===load simulation results===
load('SD_22RN_100_Mono')        %Monofacial cell
% load('SD_22RN_100_Tandem') %Tandem cell
% load('SD_22RN_100_Bifi')   %Bifacial cell
%RN     Number of rays (1*22)
%RSD    Relative Spectral Distribution (60*46)
%SEN    Sensitivity of each cell (60*22*100)
%SD     Standard deviation (100*22)
%x_co   logarithm of RN (not used)

%===module geometry used for building the data base===
l = 0.15;           %cell length [m]
w = 0.15;           %cell width [m]
L = 6;              %Unit cell length [m]
W = 1.64;           %Unit cell width [m]
RPC = l*w/(L*W)*RN; %Rays per cell

%===initialise the variables===
SEN_2  = zeros(30,size(SEN,2),size(SEN,3));
SEN_3  = zeros(20,size(SEN,2),size(SEN,3));
SEN_5  = zeros(12,size(SEN,2),size(SEN,3));
SEN_10 = zeros(6,size(SEN,2),size(SEN,3));
SEN_20 = zeros(3,size(SEN,2),size(SEN,3));
SEN_30 = zeros(2,size(SEN,2),size(SEN,3));
SEN_60 = zeros(1,size(SEN,2),size(SEN,3));

SD_2 = zeros(100,22);
SD_3 = zeros(size(SD_2));
SD_5 = zeros(size(SD_2));
SD_10 = zeros(size(SD_2));
SD_20 = zeros(size(SD_2));
SD_30 = zeros(size(SD_2));
SD_60 = zeros(size(SD_2));

%===calculate average sensitivity of each group and the standard deviation===
for u = 1:100 
    for i = 1:30    %group 2 cells
        SEN_2(i,:,u) = mean([SEN(i*2-1,:,u);SEN(i*2,:,u)]);
    end 
    SD_2(u,:) = std(SEN_2(:,:,u));
    for i = 1:20    %group 3 cells
        SEN_3(i,:,u) = mean([SEN(i*3-2,:,u);SEN(i*3-1,:,u);SEN(i*3,:,u)]);
    end
    SD_3(u,:) = std(SEN_3(:,:,u));
    for i = 1:12    %group 5 cells
        SEN_5(i,:,u) = mean(SEN((i*5-4):i*5,:,u));
    end
    SD_5(u,:) = std(SEN_5(:,:,u));
    for i = 1:6     %group 10 cells
        SEN_10(i,:,u) = mean(SEN((i*10-9):(i*10),:,u));
    end
    SD_10(u,:) = std(SEN_10(:,:,u));
    for i = 1:3     %group 20 cells
        SEN_20(i,:,u) = mean(SEN((i*20-19):(i*20),:,u));
    end
    SD_20(u,:) = std(SEN_20(:,:,u));
    for i = 1:2     %group 30 cells
        SEN_30(i,:,u) = mean(SEN((i*30-29):(i*30),:,u));
    end
    SD_30(u,:) = std(SEN_30(:,:,u));
    for i = 1       %group 60 cells
        SEN_60(i,:,u) = mean(SEN(:,:,u));
    end
    SD_60(u,:) = std(SEN_60(:,:,u));
    
end

RPC_2  = RPC*2;  RPC_3  = RPC*3;  RPC_5  = RPC*5; RPC_10 = RPC*10;
RPC_20 = RPC*20; RPC_30 = RPC*30; RPC_60 = RPC*60;

%===use all the simulation points===
SD_final = [median(SD)    median(SD_2)  median(SD_3)  median(SD_5)...
            median(SD_10) median(SD_20) median(SD_30)];
SD_final = sort(SD_final,'descend');
RPC_final = [RPC RPC_2 RPC_3 RPC_5 RPC_10 RPC_20 RPC_30];
RPC_final = sort(RPC_final,'ascend');

%===calculate number of the rays via interpolation===
rpc = spline(SD_final,RPC_final,sd);
NrRays = rpc*(L0*W0)/(l0*w0);

end