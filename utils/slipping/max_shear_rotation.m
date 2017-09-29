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
    % added rodriguez rotation (only valid for standardbasis like in cubic
    % to cubic case) myself to be independent of other matlab addons
    w = [0.   -u(3)  u(2)
         u(3)  0    -u(1)
        -u(2)  u(1)  0 ];
    w = sin(phi)*w;
    R =  (1.-cos(phi)* u(1:3)'*u(1:3) + eye(3)*cos(phi) + w   
    %vrrotvec2mat( u ); 
    %rot_originaxis_angle( rad2deg(phi), u ); --- This is only valid
    %through the origin (0,0)
    rotationVectorToMatrix(u(1:3)*phi)
    
    a = lol
    
end


end

% For symbolic computation the cross product could be written as the determinant
% of the extended basis matrix
% u = ( 1/sin(phi) ) * collect(det( cat(1,[1 1 1], m_u, ms_u) ) , g);

