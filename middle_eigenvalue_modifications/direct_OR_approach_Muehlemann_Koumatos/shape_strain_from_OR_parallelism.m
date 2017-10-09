function [ T ] = shape_strain_from_OR_parallelism( B, n, v )
% calculates transformation strains given a particular Bain strain B_i
% and an unrotated plane given by the normal vector n and a direction
% v within that plane that is unrotated.
% See Mühlemann 2016

if size(n) == [1,3]
    n = n';
    v = v';
end

Bn = cofactor( B ) * n;
m = Bn / ( sqrt( trace(Bn*Bn') );

Bv = B*v;
u = Bv / sqrt( trace(Bv*Bv');

nvc = cat(2, n, v, cross(n,v) );
muc = cat(2, m, u, cross(m,u) );

R = nvc * muc;

T = R*B;

end

