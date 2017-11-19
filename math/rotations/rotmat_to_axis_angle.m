function [ angle, axis ] = rotmat_to_axis_angle( R )
% call: [ angle, axis ] = axis_angle_to_rotmat( R )
% return the axis angle representation ([1x3 unit vec], angle in degree) of a general 3D rotation

% moved this check to function - signed_angle_from_rotmatrix
% if abs(det(R)- 1.) < 1.e-5
%     angle = acos( (trace(R) - 1.) / 2. );
% else
%     error('input is neither Identity nor a rotation matrix - please fix...')
% end

% PET:17.11.17
angle = signed_angle_from_rotmatrix( R );

axis = 1./(2*sin(angle));
axis = axis * [R(3,2)-R(2,3)  ,  R(1,3)-R(3,1)  ,  R(2,1)-R(1,2)] ;

% covert angle to degree
angle = rad2deg(angle);

end

