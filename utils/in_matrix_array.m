function bool = in_matrix_array( mat, mat_array)
% this function takes a NxN matrix - mat and an NxNxP matrix array - mat_array
% and checks whether the matrix mat is contained in mat_arrays

bool = false;
for j = 1:size(mat_array,3)
    if abs(mat - mat_array(:,:,j) ) < 1.e-5
        bool = true;
        break
    end
end

end

