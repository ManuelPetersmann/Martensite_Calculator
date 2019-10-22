function [ OR ] = OR_from_shape_strain( ST, e_i )
% given a shape_strain ST and the Bain axis of the Transformation e_i (e.g.
% [0 0 1] ) calculates the Orientation relation ship matrix (rotation)

CP = axis_angle_to_rotmat(-45. , e_i ); %rot_originaxis_angle( -45. , e_i );

[~,R] = polardecomposition( ST );

OR = R * CP;  % active OR
% the inverse gives the passive OR = coordinate transformation matrix
OR = CP' * R'; 

end

