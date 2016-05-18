function [ S ] = mirror_by_plane( n, basis1 )
% basis1 - basis vectors as colums
% given a plane normal n computes the mirror matrix for this plane.
% Note: the colums of a mapping are the transformations of the respective basis vectors

%S = [cos(2*a), sin(2*a); 
%     sin(2*a), cos(2*a)]

% h holds the scalar products of the basis vektors and the plane normal
h = basis1'*n;
n2 = dot(n,n);
lambda = [0, 0, 0];
S = zeros(3);

for i = 1:3
        lambda(i) = -h(i)/n2;
        S(:,i) = basis1(:,i) + 2*lambda(i)*n;
end

end
