function reduced = reduce_vecs( in )
% this function takes an NxP matrix of 1xN vectors and sorts out 
% entries occuring more than once
unique_entries = 0;
reduced = [0. 0. 0.];
for i = 1:size(in,1)
    if ~in_vec_array( in(i,:), reduced )
        unique_entries = unique_entries +1;
        reduced(unique_entries,:) = in(i,:);
    end
end

end

