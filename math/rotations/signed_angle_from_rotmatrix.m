function phi = signed_angle_from_rotmatrix( R )
% call: angle_from_rotmatrix(R)
% given a rotation matrix this function returns its signed rotation angle
% in radians

if abs(det(R)- 1.) > 1.e-5
    error('input is neither Identity nor a rotation matrix - please fix...')
end

abs_phi = acos( (trace(R)-1.) / 2.);
n = rand(3,1); % Note n must not be equal to axis of rotation

% rotation axis "u" is the eigenvector of R with corresponding eigenvalue 1.
[ y1, y2, y3, e1, e2, e3] = sorted_eig_vals_and_vecs( R, false ); % false -> to get no error for complex eigenvalues
yyy = [y1,y2,y3];
eee = [e1';e2';e3'];
if any(abs(yyy-1.) < 1.e-6)
    idx = find(abs(yyy-1.) < 1.e-6);
    eee(idx,:)
    sig = sign( dot( cross(n, R*n), eee(idx,:) ) ) ;
    phi = sig * abs_phi;
else
    error('Rotation matrix has no eigenvalue of 1 are you sure it is right?')
end


end

