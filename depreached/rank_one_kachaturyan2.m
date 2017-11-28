function [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one_kachaturyan2( F )
% solves the equation for an invariant habit plane to an undeformed state:
% Q*A - I = eps_0* a \otimes n 
% where B is the (symmetric) Bain strain, R is a rotation
% "a" is the normalized shape deformation direction, times length eps_0
% "h" is the normal to the aust-mart habit plane
% see Bhattacharya - Qi, Khachaturyan 2014 - Dislocated microstructure of
% martensite steel - Acta Materilia

[ y1, y2, y3, e1, ~, e3 ] = sorted_eig_vals_and_vecs( F'*F );
%'should be normalized eigenvectors KH, norm(e1), norm(e3) =...'
%norm(e1)
%norm(e3)

if (y1 > 1.  ||  abs(1-y2) > 1.e-2  ||  y3 < 1.)
    error('Eigenvalues do not satisfy conditions y1<1 , y2=1 , y3>1 necessary for an invariant plane')
end

eps_0 = sqrt(y3) - sqrt(y1);

ah1 = sqrt( y3 * (1- y1) / (y3-y1) ) * e1;  
ah3 = sqrt( y1 * (y3 -1) / (y3-y1) ) * e3; 
% 
hh1 = sqrt( (1-y1) / (y3-y1) )*e1;
hh3 = sqrt( (y3-1) / (y3-y1) )*e3;
%
% there are two solutions (different sign!)
h1 = hh3 + hh1;
h2 = hh3 - hh1;
%
a1 = ah3 - ah1;
a2 = ah3 + ah1;
%
Q1 = (eye(3) + eps_0*(a1 * h1'))*inverse(F);
Q2 = (eye(3) + eps_0*(a2 * h2'))*inverse(F);
%
end


