function bool = in_array( subarray, array, dim, negatives)
% call: in_array(subarray, array, dim, negatives)
% this function takes an array "subarray" and searches in the given dimension
% "dim" if this subarray is in the "array" e.g. a 1x3 vector in a Nx3 matrix
% "negatives" is a boolean specifying if also subarrays multiplied by (-1)
% should be sorted out (true) e.g. because the vector [1 2 3] equals [-1 -2 -3]
% programmed for array of matrices and vectors

if nargin < 4
    negatives = true;
end

bool = false;
for j = 1:size(array, dim)
    % determine size of subarray for comparison
    if size(subarray) == [1,3]
            entry = array(j,:);
    end
    if size(subarray) == [3,3]
            entry = array(:,:,j);
    end         
    %        
    if negatives == true
        if ( sum(abs(subarray - entry )) < 1.e-4 ) || ( sum(abs(subarray + entry )) < 1.e-4 )
            bool = true;
            break
        end
    elseif abs(subarray - entry ) < 1.e-2
        bool = true;
        break
    end
end

end

