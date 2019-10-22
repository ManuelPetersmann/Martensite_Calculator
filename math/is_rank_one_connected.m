function [bool,y2]  = is_rank_one_connected(F,G,tolerance)

bool = true;

Ct = inverse(G)'*(F'*F)*inverse(G); 

% no habit plane in this case
if abs(Ct - eye(3)) < tolerance
    bool = false;
    return
end

% otherwise automatically the eigenvalues are all positive!
[ y1, y2, y3] = sorted_eig_vals_and_vecs( Ct );

if (y1 > 1. ||  abs(1-y2) > tolerance  ||  y3 < 1.)
    bool = false;
    return
end


end

