function bool = linear_dependent( B1 )
% given the basisvectors as colums of a matrix B1, calculates
% wheter the basis is linear independent
if abs( det(B1) ) < 1.e-9
    bool = true;
else
    bool = false;
end
end