function cof = cofactor( A )
% Calculates the cofactor of matrix A
% surface elements transform according to cof(F)
detA = det(A);
if abs( detA ) > 1.e-9 % unequal zero    
    cof = detA*inverse(A)' ; %det(A)*A^{-T}
else
    error('matrix not invertible det uneq 0')
end
end

