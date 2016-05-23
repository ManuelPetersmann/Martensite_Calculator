
% highly symmetric mirror planes from bcc
% {001} family
m(:,1) = [0. 0. 1.];
m(:,2) = [0. 1. 0.];
m(:,3) = [1. 0. 0.];
% {011} family
m(:,4) = [0. 1. 1.];
m(:,5) = [1. 0. 1.];
m(:,6) = [-1. 0. 1.];
m(:,7) = [1. 1. 0.];
m(:,8) = [0. -1. 1.];
m(:,1) = [1. -1. 0.];

for im=1:size(m,2) % number of considered mirror planes in martensite
    
    m_mart = m(:,im); % mirror plane in martensite
    m_aust = inverse(cp) * m_mart;
    
    for is = 1:size(d,2) %Anzahl slip systems
        
        d1 = d(:,is);
        n1 = n(:,is);
        d2 = mirror_by_plane(d1,eye(3));
        n2 = mirror_by_plane(n1,eye(3));
        
        % S1 and S2 are shears related by mirror symmetry
        S1 = get_shear_d_n( d1, n1, m);
        S2 = get_shear_d_n( d2, n2, m);
        
        m_aust_sheared = inverse(S)* m_mart;
        
        R = max_shear_rotation(m_aust, S);
        
        A = 0.5*( inverse(R)*S1 + R*S2 )* B3; % composite block deformations
        
    end
    
    
end

end