function [epsilon, a1, a2, n1, n2, Q1, Q2] = rank_one(A, B)
% solves the equation for an invariant interface:
% QA - B = a \otimes b
% where A and B are two mappings (homogeneous deformations)
% "a" is the shape deformation direction
% "n" is the normal to the invariant interface
% see Bhattacharya - Microstructure of Martensites p.69

% check if this gives the same result as that of khachaturyan 
% if it does then this is more general since here the vectors
% are normed automatically

Ct = inv(B)'*(A'*A)*inv(B); 
if abs(Ct - eye(3)) < 1.e-4
    error('There is no solution')
end
% otherwise automatically the eigenvalues are all positive!

[V,D] = eig(A);
% D... Diagonal matrix of eigenvalues
% V... colum vectors are correspondig right eigenvectors i.e. A*V = V*D
eigs = [D(1,1),D(2,2),D(3,3)]
[eigs_sort, old_idx] = sort(eigs)
% arrange eigenvectors to order of sorted eigenvectors
y1 = eigs_sort(1)
y2 = eigs_sort(2)
y3 = eigs_sort(3)
epsilon = sqrt(y3) - sqrt(y1);
e1 = V(:,old_idx(1));
%e2 = V(:,old_idx(2));
e3 = V(:,old_idx(3));

if y1 > 1  ||  abs(1-y2) > 1.e-12
    error('Eigenvalues do not satisfy conditions y1<1 , y2=1 , y3>1 necessary for an invariant plane')
end

a1 = sqrt((y3*(1-y1))/(y3-y1))*e1; 
a3 = sqrt(y1*(y3-1)/(y3-y1))*e3;
% 
nh1 = (-sqrt(1-y1))*e1*B';
nh3 = X*sqrt(y3-1)*e3*B';
%
coef = ((sqrt(y3)-sqrt(y1)) / (sqrt(y3-y1)));
% there are two solutions (different sign!)
nn1 = coef* (nh1 + nh3);
nn2 = coef* (nh1 - nh3);
%
% rho unequal 0 is a constant such that |n| =1 
roh1 = norm(nn1,2);
roh2 = norm(nn2,2); % euclidic Norm
%
a1 = roh1 * (a1 + a3);
a2 = roh2 * (a1 - a3);
%
n1 = nn1/roh1;
n2 = nn2/roh2;
%
% Rotations relationg the deformations on either side
% If B = Identity then determines how a thin plate of the crystal created
% the invariant line strain "B*R" or "F" or "A0" should be reoriented to 
% achieve an invariant planar match with the parent.
Q1 = ((a1 * n1') + B)*inv(A);
Q2 = ((a2 * n2') + B)*inv(A);
%
end

