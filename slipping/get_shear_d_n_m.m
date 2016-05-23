function shear = get_shear_d_n( d, n, m )
% given the slip direction d and the slip plane normal n and the shear magnitude
% 1/m (i.e. the lattice has a step after each m planes after the shear)
% get the shear matrix in the same base

shear = eye(3) + (1./m)*(d * n'); % last term is kronecker product

end

