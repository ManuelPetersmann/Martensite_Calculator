function h_res = hp_two_IPS_same_plane(h, d1, d2, eps1, eps2)
% see Bhadeshia - worked examples in the geometry of crystals page 49 +
% formelsammlung tensoralgebra particularly dyads - on WIKI

h = h / norm(h);
d2 = d2 / norm(d2);
d1 = d1 / norm(d1);

r = eps1*d1 + eps2*d2 + eps1*eps2*dot(h*d2)*d; 

h_res = eye(3) + kron( r , h );

end

