function [index_matrix] = Find_empyty(Str_cell)
index_matrix=zeros(size(Str_cell));
for xsite=1:size(Str_cell,1)
    for ysite=1:size(Str_cell,2)
        if isempty(Str_cell{xsite,ysite})
            index_matrix(xsite,ysite)=0;
        else
            index_matrix(xsite,ysite)=1;
        end
    end
end
end