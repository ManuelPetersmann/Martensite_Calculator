function R = max_shear_rotation( m, S )
% given the shear plane normal n and the shear matrix S, this function
% calculates the rotation matrix for the vector n due to the shear.

ms = S*m;
phi = acos( dot(m,ms) / (norm(m)*norm(ms)) ); % verify notation for symbolic

m_u = m / norm(m);
ms_u = ms / norm(ms);

% write the cross product as the determinant of the extended basis matrix
% u = ( 1/sin(phi) ) * collect(det( cat(1,[1 1 1], m_u, ms_u) ) , g);  
u = ( 1./sin(phi) ) * cross( m_u , ms_u );

R = rot_originaxis_angle( phi, u );

end

