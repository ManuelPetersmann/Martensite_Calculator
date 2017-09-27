function R = max_shear_rotation( m, S )
% call: R = max_shear_rotation( m, S )
% given the shear plane normal m (colum) and the shear matrix S, this function
% calculates the rotation matrix for the vector m due to the shear.

ms = inverse(S)' * m;

phi = abs( acos( dot(m,ms) / (norm(m)*norm(ms)) ) );
% rad2deg(phi)

if rad2deg(phi) > 90 % cos(phi) < 0
    error('A simple shear cannot produce an angle >90, error! check code!')
end

if phi < 1.e-8
    R = eye(3);
else
    m_u = m / norm(m);
    ms_u = ms / norm(ms);
    u = ( 1./sin(phi) ) * cross( ms_u , m_u ); % sense of rotation is from ms to m
    % Rechtssystem via cross product ensured
    u(4) = phi;
    R =  vrrotvec2mat( u ); %rot_originaxis_angle( rad2deg(phi), u );
end


end

% For symbolic computation the cross product could be written as the determinant
% of the extended basis matrix
% u = ( 1/sin(phi) ) * collect(det( cat(1,[1 1 1], m_u, ms_u) ) , g);

