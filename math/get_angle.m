function [ theta ] = get_angle( vec1, vec2 )
% call: get_angle( vec1, vec2)
% calculation of angle theta (in degree) between two vectors in 3D space 
% (e.g. normal vectors of two planes, in order to determine misorientation angle)
% vrrotvec - gives rotation matrix between vectors

theta = abs( acosd( dot (vec1,vec2) / (norm(vec1)*norm(vec2)) ) ); %acosd - degree
if theta > 90
    theta = 180. - theta;
end

end

