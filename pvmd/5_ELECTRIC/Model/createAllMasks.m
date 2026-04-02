function [MASK,cfg] = createAllMasks(x,nc_r)
%This function creates masks for a 96-cell module. The exact layout of the
%module is shown below - each number signifies 1 cell on the module.
%The value of each digit shows which cell block does the cell belongs to.
%  | 1 1 1 1 | 2 2 2 2 | 3 3 3 3 |
%  | 1 1 1 1 | 2 2 2 2 | 3 3 3 3 |
%  | 1 1 1 1 | 2 2 2 2 | 3 3 3 3 |
%  | 1 1 1 1 | 2 2 2 2 | 3 3 3 3 |
%   -----------------------------
%  | 4 4 4 4 | 5 5 5 5 | 6 6 6 6 |
%  | 4 4 4 4 | 5 5 5 5 | 6 6 6 6 |
%  | 4 4 4 4 | 5 5 5 5 | 6 6 6 6 |
%  | 4 4 4 4 | 5 5 5 5 | 6 6 6 6 |
%
%Output: masks for the 27 possible configurations for this module
%Author - Andres Calcabrini (added to PVMD toolbox by Devyani Salokhe)
    MASK(:,:,1) =  [1*x 2*x 3*x; 4*x 5*x 6*x]; %s1p6
                
    MASK(:,:,2) =  [1*x 1*x 2*x; 2*x 3*x 3*x]; %s2p3
                
    MASK(:,:,3) =  [1*x 1*x 2*x; 3*x 2*x 3*x]; %s2p3
                
    MASK(:,:,4) =  [1*x 1*x 2*x; 3*x 3*x 2*x]; %s2p3
                
    MASK(:,:,5) =  [1*x 2*x 1*x; 2*x 3*x 3*x]; %s2p3
                
    MASK(:,:,6) =  [1*x 2*x 1*x; 3*x 2*x 3*x]; %s2p3
                
    MASK(:,:,7) =  [1*x 2*x 1*x; 3*x 3*x 2*x]; %s2p3
                
    MASK(:,:,8) =  [1*x 2*x 2*x; 1*x 3*x 3*x]; %s2p3
                
    MASK(:,:,9) =  [1*x 2*x 3*x; 1*x 2*x 3*x]; %s2p3
                
    MASK(:,:,10) = [1*x 2*x 3*x; 1*x 3*x 2*x]; %s2p3
                
    MASK(:,:,11) = [1*x 2*x 2*x; 3*x 1*x 3*x]; %s2p3
                
    MASK(:,:,12) = [1*x 2*x 3*x; 2*x 1*x 3*x]; %s2p3
                
    MASK(:,:,13) = [1*x 2*x 3*x; 3*x 1*x 2*x]; %s2p3
                
    MASK(:,:,14) = [1*x 2*x 2*x; 3*x 3*x 1*x]; %s2p3
                
    MASK(:,:,15) = [1*x 2*x 3*x; 2*x 3*x 1*x]; %s2p3
                
    MASK(:,:,16) = [1*x 2*x 3*x; 3*x 2*x 1*x]; %s2p3
                
    MASK(:,:,17) = [1*x 1*x 1*x; 2*x 2*x 2*x]; %s2p3
                
    MASK(:,:,18) = [1*x 1*x 2*x; 1*x 2*x 2*x]; %s3p2
                
    MASK(:,:,19) = [1*x 1*x 2*x; 2*x 1*x 2*x]; %s3p2
                
    MASK(:,:,20) = [1*x 1*x 2*x; 2*x 2*x 1*x]; %s3p2
                
    MASK(:,:,21) = [1*x 2*x 1*x; 1*x 2*x 2*x]; %s3p2
                
    MASK(:,:,22) = [1*x 2*x 1*x; 2*x 1*x 2*x]; %s3p2
                
    MASK(:,:,23) = [1*x 2*x 1*x; 2*x 2*x 1*x]; %s3p2
                
    MASK(:,:,24) = [1*x 2*x 2*x; 1*x 1*x 2*x]; %s3p2
                
    MASK(:,:,25) = [1*x 2*x 2*x; 1*x 2*x 1*x]; %s3p2
                
    MASK(:,:,26) = [1*x 2*x 2*x; 2*x 1*x 1*x]; %s3p2
                
    MASK(:,:,27) = [1*x 1*x 1*x; 1*x 1*x 1*x]; %s6p1
    
                
    %cfg contains 27 matrices
    %the numbers correspond to the different blocks in the module
    %the module is numbered like this
    %-------------------
    %|  1  |  2  |  3  |
    %|  4  |  5  |  6  |
    %-------------------
    %the numbers in the same row in cfg correspond to blocks connected in
    %series.
    
    cfg{1} = [1; 2; 3; 4; 5; 6];
    cfg{2} = [1 2; 3 4; 5 6];
    cfg{3} = [1 2; 3 5; 4 6];
    cfg{4} = [1 2; 3 6; 4 5];
    cfg{5} = [1 3; 2 4; 5 6];
    cfg{6} = [1 3; 2 5; 4 6];
    cfg{7} = [1 3; 2 6; 4 5];
    cfg{8} = [1 4; 2 3; 5 6];
    cfg{9} = [1 4; 2 5; 3 6];
    cfg{10} = [1 4; 2 6; 3 5];
    cfg{11} = [1 5; 2 3; 4 6];
    cfg{12} = [1 5; 2 4; 3 6];
    cfg{13} = [1 5; 2 6; 3 4];
    cfg{14} = [1 6; 2 3;4 5];
    cfg{15} = [1 6; 2 4; 3 5];
    cfg{16} = [1 6; 2 5; 3 4];
    cfg{17} = [1 2 3; 4 5 6];
    cfg{18} = [1 2 4; 3 5 6];
    cfg{19} = [1 2 5; 3 4 6];
    cfg{20} = [1 2 6; 3 4 5];
    cfg{21} = [1 3 4; 2 5 6];
    cfg{22} = [1 3 5; 2 4 6];
    cfg{23} = [1 3 6; 2 4 5];
    cfg{24} = [1 4 5; 2 3 6];
    cfg{25} = [1 4 6; 2 3 5];
    cfg{26} = [1 5 6; 2 3 4];
    cfg{27} = [1 2 3 4 5 6];
    
    if nc_r ~= 8
       MASK = permute(MASK,[2 1 3]);
%        cfg = cellfun(@transpose,cfg,'UniformOutput',false);
    end
end