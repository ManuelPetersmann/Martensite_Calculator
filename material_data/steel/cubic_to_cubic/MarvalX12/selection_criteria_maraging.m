
% densest packed plane in austenite
cpps_gamma = all_from_family_perms( [1 1 1] );
austenite.CPPs = cpps_gamma;

% % 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
% densest packed direction in austenite
KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
austenite.CP_dirs = KS;

% Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [-110]_aust || [100]_mart';
NW = all_from_family_perms( [1 2 1], false );

%% lath level selection criteria
% numerical lower bound is 5!
g_min = 5.; %  % could also directly be specified in mod_eigenvalue function e.g. block_symmetric_shear

% NOTE also check absolute bounds of g_min and eps_max in numerical
% parameters!!!!!!

% maximum value of shape strain of deformation (lambda1-lambda3) - property
% name is eps_ips
% check also set numerical upper bound 
eps_max = 0.6; %1 0.3; % 100.; 

% I put this value to 2 considering that in bcc there is no close packed  plane
theta_CPPs_max = 5.; %3.; %1.; % 90.;  maximum misorientation angle of CP relation - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice

% maximum misorientation angle of block habit-plane to {111}_gamma
theta_h_to_CPP = 20.; %10. % normally between 10 and 20 - see Maresca paper.;  
% on the block level this tolerance is reduced!
% note 557 is 9.4 degree from 111 ! therefore this high tolerance!

% The peak of OR distribution is normally between KS and NW and
% these two are 5.25 apart - hence these tolerances
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )
theta_KS_max = 6.; % 5.; %3.5; % 10.;  had 6 here in last calculations

%'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
theta_NW_max = 8.; %5.5; %3.5; % what is allowed, theta_NW_min 
% what is the smallest angle within the family-transformed_family set

delta_determinant_max = 0.001;

% PET 10.10.17
% angle between preferred ILS direction and habit plane (0 if vector in
% habit plane)
theta_max_ILSdir_to_h = 3.; % 3 reduces it to only 16 from 872 (last reduction step) with criteria from 27.10.17
%theta_e1_cpp_dir = 5.; angle between direction with smallest deformation (e1) and
% close packed direction in austenite NOT REASONABLE!!!

% Angle between 111 and 557 habit plane
% acos( dot([1. 1. 1.], [5. 5. 7.])/(sqrt(3)*sqrt(99) ) ) = 9.4 degree

%display(['Selection criteria: g_min a']) 

%% Block level selection criteria
%theta_hps = 10;
%theta_intersec_cpdir = 10.;

% block tolerances
rot_angle_block = 3.
lambda2_tol_block_aust = 1.e-3 % doesnt matter if 0.001 or 0.0001 !!! important! some more solutions with 0.003
block_hp_cp_aust_tol = 5.; % degree - even if i just set this only to 10 most solutions fall out
%lambda2_tol_laths = 1.e-4

%% Packet level selection criteria
% it seems reasonable that each packet has approximately a hydrostatic 