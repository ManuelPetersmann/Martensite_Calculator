clear all; clc;

a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 

Bain_and_Correspondence;

count_directions_extra = true;

%% assemble slip systems in alpha
% since the shear is a substantial part of the transformation only 
% shear systems which are favorable in the b.c.c. lattices are considered. 
% the plane and direction families are {110}_alpha, {112}_alpha,
% <111>_alpha, <110>_alpha
plane_families_bcc =     [ [1 1 0]
                           [1 1 2] ];   % must be written with linebreak or ";" between vectors!                     
direction_families_bcc = [ [1 1 1]
                           [1 1 0] ];
% find all possible combination (including different shear directions)
[martensite.slip_planes, martensite.slip_directions] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );
%[ ns_product, ds_product ] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );

plane_families_fcc =     [ [1 1 1] ];
direction_families_fcc = [ [1 1 0]; ] %; [1 1 2] ]; % 112 is twinning dislocation... - 112 is a partial! do not take it!
[austenite.slip_planes, austenite.slip_directions] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);
%[ ns_parent, ds_parent] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);

martensite.considered_plasticity = 3; % both mart and aust slip systems
%% calculate possible solutions and store solution objects in an object array
all_sols = multiple_shears_incremental_optimization( martensite, austenite);


%% further checks if solution is appropriate - reduction of total solutions one at a time
% criteria for selection of solutions:
% -) angular deviation to nearest cpp theta_n = min{ angle(n, {111} ) } - must be small
% -) slip density parameter g... number of planes between - steps - must be high
% -) OR's defiations: theta_p = min[ {111}_gamma < {011}_alpha = AL^-T * {111}_gamma } - criterion NR_1 see Qi2013 p.28
%    theta_KS and theta_NW
% -) shape strain - eps_0 - must be small
% -) determinant must be invariant (could change due to additive mixture of matrices

% all selection parameters have been moved to a seperate file to enable
% easier comparison between methods - here they have been commented.
% load all parameters from the file:
%selection_criteria_maraging;

%% reduce solutions to ones with g < 20. i.e. at least 20 planes between dislocations
% average number of atom layers before a step due to the (continuum) applied shear occurs (LIS)
%g_min = 10.; % 5.; % could also directly be specified in mod_eigenvalue function e.g. block_symmetric_shear
%% g_min_sols = Solution_array( Slip_solution_doubleshear(), theta_p_sols, 'g', g_min, 'min'); 


%% reduce soltuions to ones with eps < something 
%eps_max = 0.3; % 100.; 
%% eps_max_solutions = Solution_array( Slip_solution_doubleshear(), g_min_sols, 'eps', eps_max, 'max' ); 


%%
%cpps_gamma = all_from_family_perms( [1 1 1] );
% 'misorientation of c.p.p martensite to austenite';
%theta_p_max = 5.; % 90.;  maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice
%% theta_p_sols = Solution_array( Slip_solution_doubleshear(), det_sols, cpps_gamma, theta_p_max, 'theta_p', 'closest_to_cpp', 'cpps_gamma', true);


%% 
%theta_n_max = 2.; % 90.; % % maximum misorientation angle of block habit-plane to {111}_gamma
% specify family near to which habit plane solutions should be searched
% calculation of theta_n - deviation of solution from {111}
%% cpp_deviation_sols = Solution_array( Slip_solution_doubleshear(), eps_max_solutions, cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h'); 
% cpp_deviation_sols = Solution_array( Slip_solution_doubleshear(), all_sols, cpps_gamma, theta_n_max, 'theta_a', 'closest_to_a', 'a'); 
% alternatively {557}_gamma could be used here see Iwashita 2011


%% specify for directions the misorientation should be small in the solutions
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
%theta_KS_max = 3.; % had 6 here in last calculations
%KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )
%% theta_KS_sols = Solution_array( Slip_solution_doubleshear(), cpp_deviation_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );


% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
%% theta_NW_sols = Solution_array( Slip_solution_doubleshear(), theta_KS_sols, NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);

%% det_sols = Solution_array( Slip_solution_doubleshear(), all_sols, 'det', delta_determinant_max,  detB3); 
%
% to sort fully reduced solution for most important criterion 
%% qi_sols = det_sols.sort( 'theta_p' ); % sort in ascending order for specific property
%theta_NW_sols.array(1) % print out best solution



  



       



