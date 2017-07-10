function vec_mirror = mirror_by_plane( n, vec, basis1 )
% call: mirror_by_plane( n, vec, basis)
% "basis1" - basis vectors as colums
% given a plane_normal "n" (colum vector), first computes the mirror matrix for this plane.
% and then the mirror image of a vector "vec" (colum vector).

% Note: the colums of a mapping are the transformations of the respective basis vectors

% h holds the scalar products of the basis vectors and the plane normal
% (Basiszerlegung)
h = basis1'*n;

lambda = [0, 0, 0];
S = zeros(3);
for i = 1:3
        lambda(i) = -h(i)/dot(n,n);
        S(:,i) = basis1(:,i) + 2*lambda(i)*n;
end

vec_mirror = S*vec;

end
