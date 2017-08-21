cpps_gamma = all_from_family_perms( [1 1 1] );
KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
NW = all_from_family_perms( [1 2 1], false );

%g_min = 10.; % 5.; % could also directly be specified in mod_eigenvalue function e.g. block_symmetric_shear

eps_max = 0.3; % 100.; 

theta_p_max = 1.; % 90.;  maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice

theta_n_max = 10.; % 90.; % % maximum misorientation angle of block habit-plane to {111}_gamma
% note 557 is 9.4� from 111 ! therefore this high tolerance!


% The peak of OR distribution is normally between KS and NW and
% these two are 5.25 apart - hence these tolerances
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
theta_KS_max = 3.5; % had 6 here in last calculations
%'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
theta_NW_max = 3.5; % what is allowed, theta_NW_min - what is the smallest angle within the family-transformed_family set

delta_determinant_max = 0.0001;

% Angle between 111 and 557 habit plane
% acos( dot([1. 1. 1.], [5. 5. 7.])/(sqrt(3)*sqrt(99) ) ) = 9.4 degree
