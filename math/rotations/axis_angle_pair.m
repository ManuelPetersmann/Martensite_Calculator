function [ angle, axis ] = axis_angle_pair( R )
% ANGLE_FROM_ROTMATRIX
% return the angle of a general 3D rotation in degrees

if abs(det(R)- 1.) < 1.e-4
    angle = acos( (trace(R) - 1.) / 2. );
else
    error('input is not a rotation matrix - please fix...')
end

axis = 1./(2*sin(angle));
axis = axis * [R(3,2)-R(2,3)  ,  R(1,3)-R(3,1)  ,  R(2,1)-R(1,2)] ;

% covert angle to degree
angle = rad2deg(angle);

end

