function phi = signed_angle_from_rotmatrix( R )
% call: angle_from_rotmatrix(R)
% given a rotation matrix this function returns its signed rotation angle

abs_phi = acos( (trace(R)-1.) / 2.);
n = rand(3,1); % Note n must not be equal to axis of rotation

% rotation axis "u" is the eigenvector of R with corresponding eigenvalue 1.
[ ~, y2, ~, ~, u, ~] = sorted_eig_vals_and_vecs( R );

if (y2 - 1.) < 1.e-6
    sig = sign( dot( cross(n, R*n), u ) ) ;
    phi = sig * abs_phi;
else
    error('Rotation matrix has no eigenvalue of 1 are you sure it is right?')
end


end

