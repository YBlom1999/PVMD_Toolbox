function S = concatstruct(S1,S2)
%combine 2 structs which may have different sets of fields

ff = unique([fieldnames(S1);fieldnames(S2)]);   %get ALL fieldnames
for n = 1:length(ff)                            %for every name
    if ~isfield(S1,ff(n)), S1(1).(ff{n})=[]; end   %if need add to S1
    if ~isfield(S2,ff(n)), S2(1).(ff{n})=[]; end   %if need add to S2
end
S = [S1,S2];  %standard concat requires same set of field names
end