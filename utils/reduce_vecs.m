function reduced = reduce_vecs( vecs, negatives )
% call: reduce_vecs( vecs )
% this function takes an Nx3 matrix of 1x3 vectors and sorts out 
% vectors occuring twice, i.e. vec1 = vec2 and vec1 = (-1)*vec2
% "negatives" is a boolean specifying if also subarrays multiplied by (-1)
% should be sorted out (true) e.g. because the vector [1 2 3] equals [-1 -2 -3]

if nargin < 2
    negatives = false; % edited to false ---  29.08.17 Manuel
end

unique_entries = 0 ;
reduced = [0. 0. 0.];

for i = 1:size(vecs,1)
    if ~in_array( vecs(i,:), reduced, 1, negatives)
        unique_entries = unique_entries +1;
        reduced(unique_entries,:) = vecs(i,:);
    end
end

end

