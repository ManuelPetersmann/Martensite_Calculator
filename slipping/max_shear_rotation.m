function R = max_shear_rotation( n, S )
% given the shear plane normal n and the shear matrix S, this function
% calculates how much the vector n rotates due to the shear.

ns = S*n;
phi = acos( dot(n,ns) / (norm2(n)*norm2(ns)) );

n_u = n / norm2(n);
ns_u = ns / norm2(ns);

u = ( 1/sin(phi) ) * cross( n_u , ns_u );

R = rot_originaxis_angle( phi, u )

end

