function C = contraction_2_1( A, B )
% call: contraction_2_1( B, C )
% contracts two second order tensors B and C to a scalar

C = 0;
for i = 1:3
    for j = 1:3
        C = C + A(i,j) * B(i,j);
    end
end

end

