function [ angle ] = angle_from_rotmatrix( R )
% ANGLE_FROM_ROTMATRIX
% return the angle of a general 3D rotation in degrees

if abs(det(R)- 1.) < 1.e-4
    angle = acosd( (trace(R) - 1.) / 2. );
else
    error('input is not a rotation matrix - please fix...')
end

end

