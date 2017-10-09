% densest packed plane in austenite
cpps_gamma = all_from_family_perms( [1 1 1] );
% % 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
% densest packed direction in austenite
KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
% Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [-110]_aust || [100]_mart';
NW = all_from_family_perms( [1 2 1], false );

%% lath level selection criteria
g_min = 5.; %  % could also directly be specified in mod_eigenvalue function e.g. block_symmetric_shear

% maximum value of shape strain of deformation (lambda1-lambda3) - property
% name is eps_ips
eps_max = 0.6; %1 0.3; % 100.; 

theta_CPP_max = 1.5; %1.; % 90.;  maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice

% on the block level this tolerance is reduced!
theta_n_max = 11.; % normally between 10 and 20 - see Maresca paper.; % 
% maximum misorientation angle of block habit-plane to {111}_gamma
% note 557 is 9.4° from 111 ! therefore this high tolerance!

% The peak of OR distribution is normally between KS and NW and
% these two are 5.25 apart - hence these tolerances
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )
theta_KS_max = 3.5; % 10.;  had 6 here in last calculations
%'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
theta_NW_max = 5.5; %3.5; % what is allowed, theta_NW_min 
% what is the smallest angle within the family-transformed_family set

delta_determinant_max = 0.0001;

% Angle between 111 and 557 habit plane
% acos( dot([1. 1. 1.], [5. 5. 7.])/(sqrt(3)*sqrt(99) ) ) = 9.4 degree

%display(['Selection criteria: g_min a']) 

%% Block level selection criteria
theta_n_max_block = 1.;

%% Packet level selection criteria
% it seems reasonable that each packet has approximately a hydrostatic 