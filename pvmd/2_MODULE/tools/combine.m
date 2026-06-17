function [V,F,T] = combine(VV,FF,TT,BackwardTracer)
%combine vertex (V), facet (F) and type (T) of different objects
%Vertex and facet matrices can be concatenated, but this changes
%their numbering. So the facet matrix and T.Facet has to be
%renumbered, respectively. Giving TT = [] skips generation of T.

nro = length(VV);       %nr of objects

V = VV{1};              %first object part of combined vertex matrix
v = size(VV{1},1);      %vertex shift nr 1
F = FF{1};              %first object part of combined facet matrix

if ~isempty(TT)
    f = size(FF{1},1);  %face shift nr 1
    T = TT{1};              %first object part of combined type struct
else
    T = [];
end

for o = 2:nro           %for every object
    V = [V;VV{o}];      %append vertices (TODO: check for duplicates?)
    if BackwardTracer && size(FF{o},2) == 4
        F = [F;FF{o}(:,[1,2,3])+v;FF{o}(:,[1,3,4])+v];    %facet matrix points to new vertex indices
    else
        F = [F;FF{o}+v];    %facet matrix points to new vertex indices
    end
    v = v+size(VV{o},1); %increment vertex shift
    if ~isempty(TT)
        l = length(TT{o});
        for t = 1:l         %for every type
            TT{o}(t).Facet = TT{o}(t).Facet+f;
        end
        T = concatstruct(T,TT{o}); %add Type structure
        f = f+size(FF{o},1); %increment facet shift
    end
end

end