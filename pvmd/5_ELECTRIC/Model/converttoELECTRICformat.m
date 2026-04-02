
function J_new = converttoELECTRICformat(J)
J_new = cell(1,1);
for k = 1:size(J,3)
    a1 = J(:,:,k);
    a1 = rot90(rot90(a1,2));
    a1 = a1(1:end); %converting the cell values to row matrix
    J_new{1} = [J_new{1};a1];
end
end