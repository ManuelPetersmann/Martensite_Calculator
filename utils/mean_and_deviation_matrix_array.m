function [median_mat, deviation_mat] = mean_and_deviation_matrix_array( ma )
% calculates the mean matrix 1/k sum_k  m_ij of a   [u x v x k] matrix array
% and the positive and negative deviation of each entry (the larger of the
% two is given

s = size(ma);
median_mat = zeros( s(1), s(2) );
deviation_mat = zeros( s(1), s(2) );

for i = 1:s(1)
    for j = 1:s(2)
        for k = 1:s(3)
            median_mat(i,j) = median_mat(i,j) + ma(i,j,k);
        end
    end
end
median_mat = abs( median_mat / s(3) );
for i = 1:s(1)
    for j = 1:s(2)
        for k = 1:s(3)
            diff = abs( abs(ma(i,j,k)) - median_mat(i,j) );
            if diff > deviation_mat(i,j) 
                deviation_mat(i,j) = diff;
            end
        end
    end
end

end

