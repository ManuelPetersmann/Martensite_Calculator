function [ theta, closest_from_vecs] = min_misorientation( vecs, comp, plane )
% call: [theta, closest] = misorientation( vecs1, comp=vecs2/lattice_transformation, theta_max, plane)
% given an array of vectors "vecs1" representing families of directions or planes
% calculates the minimum misorientation angle "theta" (in degrees) between them and the array comp(arison) 
% and gives the family member for it. Alternatively a comp can be alattice_transformation
% matrix, used with a third argument "plane", a boolean 
% specifying wheter the mapping is applied to planes (true) or directions (false). 

theta = 999.9; % init. angle for comparison of angles between different cpp

for j = 1 : size(vecs,1)
    vec1 = vecs(j,:);
    if nargin < 3
        vec2 = comp; % comparison between family and one vector only
    elseif plane == true % nargin == 3 --> comp = LT
        vec2 = vec1 * inverse(comp);
    else
        vec2 = vec1 * comp'; % direction transformation
    end
    
    theta_new = get_angle(vec1, vec2); % returns angle in degree
    if( theta_new < theta )
        theta = theta_new;
        closest_from_vecs = vec1;
    end
end

end


