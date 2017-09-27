function [ OR ] = OR_from_shape_strain( ST, e_i )
% given a shape_strain ST and the Bain axis of the Transformation e_i (e.g.
% [0 0 1] ) calculates the Orientation relation ship matrix (rotation)

CP = rot_originaxis_angle( -45. , e_i );

[U,R] = polardecomposition( ST );

OR = CP * R'; 

end

