function [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F, G)
% solves the equation for an invariant interface:
% Q*F - G = a \otimes n
% where A and B are two mappings (homogeneous deformations)
% For an invariant interface between Austenite and one Variant of
% martensite the equation looks like this:
% Q*B - eye(3) = a \otimes n, where B is a Bain strain,
% "a" is the shape deformation direction
% "n" is the normal to the invariant interface
% see Bhattacharya - Microstructure of Martensites p.69

% check if this gives the same result as that of khachaturyan 
% if it does then this is more general since here the vectors
% are normed automatically

Ct = inverse(G)'*(F'*F)*inverse(G); 
if abs(Ct - eye(3)) < 1.e-4
    error('There is no solution')
end
% otherwise automatically the eigenvalues are all positive!
[ y1, y2, y3, e1, ~, e3 ] = sorted_eig_vals_and_vecs( Ct );

if (y1 > 1. ||  abs(1-y2) > 1.e-8  ||  y3 < 1.)
    error('Eigenvalues do not satisfy conditions y1<1 , y2=1 , y3>1 necessary for an invariant plane')
end

% shape strain eps_0
eps_0 = sqrt(y3) - sqrt(y1); 

ah1 = sqrt( y3 * (1-y1) / (y3-y1) ) * e1; 
ah3 = sqrt( y1 * (y3-1) / (y3-y1) ) * e3;
% 
hh1 =  (-1)*sqrt(1-y1) * G' *e1;
hh3 =       sqrt(y3-1) * G' *e3;
%
coef = ((sqrt(y3)-sqrt(y1)) / (sqrt(y3-y1)));
% there are two solutions (different sign!)
hh1 = coef* (hh1 + hh3);
hh2 = coef* (hh1 - hh3);%%%%%%%%%%%%
%
% rho unequal 0 is a constant such that |n| =1 
roh1 = norm(hh1,2);
roh2 = norm(hh2,2); % euclidic Norm
%
a1 = (ah1 + ah3); %roh1 * 
a2 = (ah1 - ah3); %roh2 * 
%
h1 = hh1/roh1;
h2 = hh2/roh2;
%
% Rotations relationg the deformations on either side
% If B = Identity then determines how a thin plate of the crystal created
% the invariant line strain "B*R" or "F" or "A0" should be reoriented to 
% achieve an invariant planar match with the parent.
Q1 = (G + eps_0*(a1 * h1'))*inverse(F);
Q2 = (G + eps_0*(a2 * h2'))*inverse(F);
%
end

