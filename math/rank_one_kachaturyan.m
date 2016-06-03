function [eps_0, a1, a2, n1, n2, Q1, Q2] = rank_one_kachaturyan( B )
% solves the equation for an invariant habit plane to an undeformed stat:
% I + eps_0* l \otimes n = BR = F
% where B is the (symmetric) Bain strain, R is a rotation
% "l" is the shape deformation direction
% "h" is the normal to the aust-mart habit plane
% see Bhattacharya - Qi, Khachaturyan - Dislocated microstructures
%
% Ehl: 03.06.2016: changed epsilon to eps_0 
% in order to distinguish from the epsilon as criteorion of the search-algorithm

B2 = (B'*B); 

[ y1, y2, y3, e1, e2, e3 ] = sorted_eig_vals_and_vecs( B2 );

if y1 > 1  ||  abs(1-y2) > 1.e-6
    error('Eigenvalues do not satisfy conditions y1<1 , y2=1 , y3>1 necessary for an invariant plane')
end

eps_0 = sqrt(y3) - sqrt(y1);

% Ehl: changed name here because later a1 one is calculated and then used again
ah1 = y3*sqrt( (1- y1^2)/(y3^2-y1^2) )*e1; 
ah3 = y1*sqrt( (y3^2 -1)/(y3^2-y1^2) )*e3;
% 
nh1 = sqrt( (1-y1^2)/(y3^2-y1^2) )*e1;
nh3 = sqrt( (y3^2-1)/(y3^2-y1^2) )*e3;
%
% there are two solutions (different sign!)
n1 = nh3 + nh1;
n2 = nh3 - nh1;
%
% Ehl: changed the name of variables above, because in calculation of a2, a1 was already overwritten 
% --> lead to a wrong result for a2 --> further error in calculation of Q2
a1 = ah3 - ah1;
a2 = ah3 + ah1;
%
%Ehl: rotation matrices? --> Eq.(9) ?
Q1 = (eye(3) + eps_0*(a1 * n1'))*inv(B);
Q2 = (eye(3) + eps_0*(a2 * n2'))*inv(B);
%
end


