function [H,bbox] = geotiffmerge(TIFFfiles,bbox_S)

GTindx = zeros(2,2);
for i = 1:size(TIFFfiles,1)
    bbox_tif = TIFFfiles(i).bbox;
    ix1 = bbox_tif(1,1); iy1 = bbox_tif(1,2); 
    ix2 = bbox_tif(2,1); iy2 = bbox_tif(2,2);
    for r = 1:1:2
        for c = 1:1:2
            if  ix1 <= bbox_S(r,1) && bbox_S(r,1) <= ix2 && ...
                iy1 <= bbox_S(c,2) && bbox_S(c,2) <= iy2
                    GTindx(r,3-c) = i;
            end
        end
    end
end
uGTi = unique(GTindx);
if ismember(0,GTindx)
    error('Missing LiDAR tile');
    uGTi = uGTi(uGTi~=0);
end
if numel(uGTi) == 1
    H = TIFFfiles(uGTi).H; bbox = TIFFfiles(uGTi).bbox;
else 
    disp('Merging LiDAR data...'); tic
    bbox = [min(cat(1,TIFFfiles(uGTi).bbox));max(cat(1,TIFFfiles(uGTi).bbox))];
    xlin = bbox(1,1):0.5:bbox(2,1);  % Get X vector: east->west (Low->High)
    ylin = bbox(2,2):-0.5:bbox(1,2); % Get Y vector north->south (High->Low)

    H = zeros(length(ylin),length(xlin));    
    for j=1:numel(uGTi)
        ix1 = find(xlin>=TIFFfiles(uGTi(j)).bbox(1,1),1,'first'); 
        ix2 = find(xlin<=TIFFfiles(uGTi(j)).bbox(2,1),1,'last');
        iy1 = find(ylin<=TIFFfiles(uGTi(j)).bbox(2,2),1,'first'); 
        iy2 = find(ylin>=TIFFfiles(uGTi(j)).bbox(1,2),1,'last');
        H(iy1:iy2,ix1:ix2) = TIFFfiles(uGTi(j)).H;

    disp('Merging complete!'); toc
    end          
end