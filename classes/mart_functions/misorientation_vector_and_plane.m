function [ theta, closest_from_vecs] = misorientation_vector_and_plane( vecs, plane)
% call: [theta, closest] = misorientation( vecs, plane)
% given an array of vectors "vecs" 3xN representing families of directions or planes
% and a plane normal vector 'plane' calculates the minimum misorientation angle "theta" (in degrees)
% of one of the vectors and the plane (zero if the vector lies in the plane).
% returns the angle and the vector.

theta = 999.9; % init. angle for comparison of angles between different cpp
closest_from_vecs = zeros(1,3);

for j = 1 : size(vecs,1)
    vec1 = vecs(j,:);
    theta_new = abs( 90. - get_angle(vec1, plane) ); % returns angle in degree
    if( theta_new < theta )
        theta = theta_new;
        closest_from_vecs = vec1;
    end
end

if closest_from_vecs == zeros(1,3)
    error('error in function misorientation_vector_and_plane - must be fixed!')
    % with tolerance: No solution for prefered invariant line within found habit plane with the given tolerance
end


