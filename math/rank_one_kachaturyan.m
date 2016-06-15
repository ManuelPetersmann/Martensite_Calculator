function [e0, a1, a2, n1, n2, Q1, Q2] = rank_one_kachaturyan( B )
% solves the equation for an invariant habit plane to an undeformed stat:
% I + epsilon* l \otimes n = BR = F
% where B is the (symmetric) Bain strain, R is a rotation
% "l" is the shape deformation direction
% "h" is the normal to the aust-mart habit plane
% see Bhattacharya - Qi, Khachaturyan - Dislocated microstructures

[ y1, y2, y3, e1, e2, e3 ] = sorted_eig_vals_and_vecs( B'*B );

if (y1 > 1.  ||  abs(1-y2) > 1.e-8  ||  y3 < 1.)
    error('Eigenvalues do not satisfy conditions y1<1 , y2=1 , y3>1 necessary for an invariant plane')
end

e0 = sqrt(y3) - sqrt(y1);

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
Q1 = (eye(3) + e0*(a1 * n1'))*inverse(B);
Q2 = (eye(3) + e0*(a2 * n2'))*inverse(B);
%
end


