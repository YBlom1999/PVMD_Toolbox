function tab = unslice(tab,s)
%Merges slices back together for plotting
% ---INPUT---
% s         vector indicating which columns will be summed when merging

r = 1;                  %start at row 1

while r<size(tab,1)    %r increases, table height decreases in loop

    %compare row material to next one
    if strcmp(tab.Layer(r),tab.Layer(r+1))    % if identical
        tab{r,s} = string(str2double(tab{r,s}) + str2double(tab{r+1,s}));           %sum
        tab(r+1,:) = [];                      % remove next
        %not r+1, so next iteration SAME row will be compared to NEW next
    else
        r = r + 1;         %if different material, simply move to next row
    end

end


end