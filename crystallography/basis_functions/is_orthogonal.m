function bool = is_orthogonal( obj )
% checks whether a given basis is orthogonal
bool = true;
for i = 1:(size(obj.C,2) -1)
    for j = i+1:size(obj.C,2)
        if dot(obj.C(:,i), obj.C(:,j)) > 1.e-9
            bool = false;
        end
    end
end
end

