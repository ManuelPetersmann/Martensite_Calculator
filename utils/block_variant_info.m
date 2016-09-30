function [ F1, F2, AL1, AL2 ] = block_variant_info( Slip_solution, B )
% Given a Slip_solution object - calculates the individual variants of the
% corresponding Block, i.e. their Lattice- as well as shape-transformations

I = eye(3);

d1 = Slip_solution.d1';
d2 = Slip_solution.d2';
n1 = Slip_solution.n1';
n2 = Slip_solution.n2';
m_aust = Slip_solution.m';
g = Slip_solution.g;
Q = Slip_solution.Q;


d11 = mirror_by_plane(m_aust, d1, I);
n11 = mirror_by_plane(m_aust, n1, I);
d22 = mirror_by_plane(m_aust, d2, I);
n22 = mirror_by_plane(m_aust, n2, I);

S1  = (d1  * n1') ;
S2  = (d2  * n2') ;
S11 = (d11 * n11') ;
S22 = (d22 * n22') ;

S =  I + (1./g)* (S1 + S2);
S_mirror = I + (1./g)* (S11 + S22);

% calculate the rotation of the mirror plane vector due to the shear
R = max_shear_rotation( m_aust, S);

AL1 = Q * R  * B;         % first variant
AL2 = Q * inverse(R) * B; % second variant

% Shape transformations of the two variants
F1 = Q * inverse(R)          * S * B;         % first variant
F2 = Q *        R * S_mirror * B; % second variant

end

