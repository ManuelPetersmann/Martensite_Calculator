function [epsilon, a1, a2, n1, n2, Q1, Q2] = rank_one_kachaturyan( B )
% solves the equation for an invariant habit plane to an undeformed stat:
% I + epsilon* l \otimes n = BR = F
% where B is the (symmetric) Bain strain, R is a rotation
% "l" is the shape deformation direction
% "h" is the normal to the aust-mart habit plane
% see Bhattacharya - Qi, Khachaturyan - Dislocated microstructures

B2 = (B'*B); 

[V,D] = eig(B);
% D... Diagonal matrix of eigenvalues
% V... colum vectors are correspondig right eigenvectors i.e. A*V = V*D
eigs = [D(1,1),D(2,2),D(3,3)];
[eigs_sort, old_idx] = sort(eigs);
% arrange eigenvectors to order of sorted eigenvectors
y1 = eigs_sort(1)
y2 = eigs_sort(2)
y3 = eigs_sort(3)
e1 = V(:,old_idx(1));
%e2 = V(:,old_idx(2));
e3 = V(:,old_idx(3));

if y1 > 1  ||  abs(1-y2) > 1.e-6
    error('Eigenvalues do not satisfy conditions y1<1 , y2=1 , y3>1 necessary for an invariant plane')
end

epsilon = sqrt(y3) - sqrt(y1);

a1 = y3*sqrt( (1- y1^2)/(y3^2-y1^2) )*e1; 
a3 = y1*sqrt( (y3^2 -1)/(y3^2-y1^2) )*e3;
% 
nh1 = sqrt( (1-y1^2)/(y3^2-y1^2) )*e1;
nh3 = sqrt( (y3^2-1)/(y3^2-y1^2) )*e3;
%
% there are two solutions (different sign!)
n1 = nh3 + nh1;
n2 = nh3 - nh1;
%
a1 = a3 - a1;
a2 = a3 + a1;
%
Q1 = (eye(3) + epsilon*(a1 * n1'))*inv(B);
Q2 = (eye(3) + epsilon*(a2 * n2'))*inv(B);
%
end


