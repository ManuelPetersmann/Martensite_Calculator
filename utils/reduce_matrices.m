function reduced = reduce_matrices( in )
% this function takes an NxNxP array of NxN matrices and sorts out 
% entries occuring more than once
unique_entries = 0;
reduced = zeros(3);
for i = 1:size(in,3)
    if ~in_matrix_array(in(:,:,i), reduced)
    unique_entries = unique_entries +1;
    end   
    if unique_entries == size(reduced, 3)      
        reduced(:,:,unique_entries) = in(:,:,i); 
    end
end

end

