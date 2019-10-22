function h_res = hp_two_ips_same_direction(h1, h2, d, eps1, eps2 )
% see Bhadeshia - worked examples in the geometry of crystals + page 49 +
% formelsammlung tensoralgebra particularly dyads - on WIKI

h1 = h1 / norm(h1);
h2 = h2 / norm(h2);
d = d / norm(d);

r = eps1*h1 + eps2*h2 + eps1*eps2*dot(h1*d)*h2; 

h_res = eye(3) + kron( d , r );

end

