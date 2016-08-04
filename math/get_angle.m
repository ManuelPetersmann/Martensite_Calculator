function [ theta_p ] = get_angle( vec1, vec2 )
% call: get_angle( vec1, vec2)
% calculation of angle theta (in degree) between two vectors in 3D space 
% (e.g. normal vectors of two planes, in order to determine misorientation angle)

theta_p = abs( acosd( dot (vec1,vec2) / (norm(vec1)*norm(vec2)) ) ); %acosd - degree
if theta_p > 90
    theta_p = 180. - theta_p;
end

end

