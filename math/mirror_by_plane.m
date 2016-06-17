function vec_mirror = mirror_by_plane( n, vec, basis1 )
% basis1 - basis vectors as colums
% given a plane normal n first computes the mirror matrix for this plane.
% And then calculates the mirror image of a vector vec with this matrix.
% Note: the colums of a mapping are the transformations of the respective basis vectors

% Ehl: ? where is this definition from? 
% S = [cos(2*a), sin(2*a); 
%      sin(2*a), cos(2*a)]

% h holds the scalar products of the basis vectors and the plane normal
% (Basiszerlegung)
h = basis1'*n;

n2 = dot(n,n);
lambda = [0, 0, 0];
S = zeros(3);

for i = 1:3
        lambda(i) = -h(i)/n2;
        S(:,i) = basis1(:,i) + 2*lambda(i)*n;
end

vec_mirror = S*vec;

end
