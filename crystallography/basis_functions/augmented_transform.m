
% In the usual description by columns, the
% vector coefficients cannot be distinguished from the point coordinates,
% but in the augmented-column description the difference
% becomes visible: the vector from the point P to the point Q has the
% coefficients v1 = q1?p1, v2 = q2?p2, v3 = q3?p3, 1?1. Thus,
% the column of the coefficients of a vector is not augmented by ‘1’
% but by ‘0’.
% 
% The advantage of the use of 4 x 4 matrices is that a sequence of
% affine transformations corresponds to the product of the correspond-

function basis1 = augmented_transform( P, basis1, shift )

mf = matrix_functions;
inv_P = mf.inverse( P );
M = inv_P;
M(4,4) = 1.;
M(4,1:3) = 0.;
M(1:3, 4) = (-1)*inv_P*shift; % see p.31 Ulrich Müller 

end

