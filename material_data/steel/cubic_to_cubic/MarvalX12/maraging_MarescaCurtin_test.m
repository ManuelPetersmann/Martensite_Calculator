clear all; clc;

a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 

Bain_and_Correspondence;
martensite.U = eye(3)*(1. /1.2537)*sqrt(2);
martensite.U(3,3) = martensite.U(3,3) / sqrt(2);

%% prescribed orientation relationship (angle in degree)
phi = 4.5;
count_directions_extra = true;

a_aust / a_mart
display('if above value > 1,274 then beta can be smaller one - see Maresca 2017 Acta')


%% assemble slip systems in alpha
% slip systems b.c.c  for P^(2)
% following Maresca&Curtin (2017) we consider
% the plane and direction families for P^(2) are {110}_alpha
% and <111>_alpha
plane_families_P2 =     [1 1 0];
direction_families_P2 = [1 1 1];
[ ns_P2, ds_P2 ] = independent_slipsystems(  plane_families_P2, direction_families_P2, count_directions_extra );
%[martensite.slip_planes, martensite.slip_directions] =

% slip systems b.c.c  for P^(3)
% following Maresca&Curtin (2017) we consider
% the plane and direction families for P^(3) are {111}_gamma
% and <101>_gamma
plane_families_P3 =     [ [1 1 1] ];
direction_families_P3 = [ [1 0 1] ];
[ ns_P3, ds_P3 ] = independent_slipsystems(  plane_families_P3, direction_families_P3, count_directions_extra );
% same as independent_slip_systems... 

% % % % transform from fcc to bcc
% % % ns_P3 = (inverse(cp)' * ns_P3')';
% % % ds_P3 = (cp * ds_P3')';

sort_out_negatives = false;
% highly symmetric mirror planes from bcc (not needed for calculations after Maresca&Curtin,
% but we keep them here to test the solutions of our ansatz with use of the slip-systems considered by Maresca&Curtin)
% {001} family
%ms = all_from_family_perms( [0 0 1], sort_out_negatives );
% {011} family
%ms = cat(1, ms, all_from_family_perms( [0 1 1], sort_out_negatives ) ); % This is a more reasonable plane - for a Block boundary!


% % ns = cat(1, ns_P2, ns_P3 );
% % ds = cat(1, ds_P2, ds_P3 );

%% calculate possible solutions and store solution objects in an object array
% all_sols = block_symmetric_doubleshear( B3, cp, ms, ns, ds); % for testing of Maresca&Curtins slip systems 
all_sols_MarescaCurtin = interface_defects_doubleshear_MarescaCurtin( martensite, ns_P2, ds_P2, ns_P3, ds_P3, phi); % for MARESCA&CURTIN RUN
% all_sols = block_symmetric_shear( B3, cp, ms, ns, ds);

%% further checks if solution is appropriate - reduction of total solutions one at a time
% criteria for selection of solutions:
% -) angular deviation to nearest cpp theta_n = min{ angle(n, {111} ) } - must be small
% -) slip density parameter g... number of planes between - steps - must be high
% -) OR's defiations: theta_p = min[ {111}_gamma < {011}_alpha = AL^-T * {111}_gamma } - criterion NR_1 see Qi2013 p.28
%    theta_KS and theta_NW
% -) shape strain - eps_0 - must be small
% -) determinant must be invariant (could change due to additive mixture of matrices

selection_criteria_maraging;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering according to Criterion 1: plastic accomodation
%%% [- g_min (minimal slip plane spacing)  
% Manuel 31.7.17
% here instead of g beta is saved in solution_array - calculate beta_min from g_min
% m3 = 1/g = (1/beta)*sqrt(3/2); --> beta = g*sqrt(3/2)
%
%beta_min = g_min*sqrt(1.5);
%g_min_sols = Solution_array( Slip_solution, all_sols_MarescaCurtin, 'g', beta_min, 'min'); 
%display(['criterion g_min=',num2str(g_min),'  current solutions: ' , num2str(size(g_min_sols.array,2)) ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EHL: Added: March 2017
% delta_determinant_max = 0.0001;
% det_sols = Solution_array( Slip_solution, theta_NW_sols, 'det', delta_determinant_max,  detB3 ); 
% display(['criterion delta_det_max=',num2str(delta_determinant_max),'  current solutions: ' , num2str(size(det_sols.array,2)) ]);
det_sols = Solution_array( IPS_solution, all_sols_MarescaCurtin, 'delta_determinant_max', delta_determinant_max,  det(martensite.U));
display(['with criterion tolerable volume_change_from_averaging = ',num2str(delta_determinant_max)] );

%% essential selection criteria
% Habit plane deviation from experimental observations - theta_h_to_CPP = 20.; %10. % normally between 10 and 20 - see Maresca paper.;  
%
%cpp_deviation_sols = Solution_array( Slip_solution, all_sols_MarescaCurtin, cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h'); 
%display(['criterion del_habitplane_111gamma_max=',num2str(theta_n_max),'  current solutions: ' , num2str(size(cpp_deviation_sols.array,2)) ]);
% cpp_deviation_sols = Solution_array( Slip_solution, all_sols, cpps_gamma, theta_n_max, 'theta_a', 'closest_to_d', 'd'); 
%
% alternatively {557}_gamma could be used here see Iwashita 2011
tolerable_HP_deviations = Solution_array( IPS_solution, det_sols, cpps_gamma, ...
    theta_h_to_CPP, 'theta_h_to_CPP', 'closest_to_h', 'h'); 
display(['with criterion del_habitplane_111gamma_max = ',num2str(theta_h_to_CPP)]);

%% 'misorientation of CPP martensite to austenite - planes of OR';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering according to Criterion 2: OR selection
%%% --> with help of - theta_p_max (maximum deviation of gamma- and alpha-cpps)
%%%                  - theta_KS_max (maximum deviation from KS-OR)
%%%              and - theta_NW_max (maximum deviation from NW-OR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpps_gamma = all_from_family_perms( [1 1 1] );
% 'misorientation of c.p.p martensite to austenite';
%theta_p_max = 5.; % 90.;  maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice

%theta_p_sols = Solution_array( Slip_solution, eps_max_solutions, cpps_gamma, theta_p_max, 'theta_p', 'closest_to_cpp', 'cpps_gamma', true);
%display(['criterion delta_CPP_max=',num2str(theta_p_max),'  current solutions: ' , num2str(size(theta_p_sols.array,2)) ]);
tolerable_CPP_deviations = Solution_array( IPS_solution, tolerable_HP_deviations, cpps_gamma, theta_CPPs_max, ...
    'theta_CPPs', 'closest_to_cpp', 'cpps_gamma', true);
display(['with criterion delta_CPPs_max = ',num2str(theta_CPPs_max)] );


%% specify maximum misorientations of solutions to ideal OR directions
%
% specify for directions the misorientation should be small in the solutions
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
%theta_KS_max = 10. ; %3.; % had 6 here in last calculations
%KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )

%theta_KS_sols = Solution_array( Slip_solution, theta_p_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );
%display(['criterion KS_direction_deviation_max=',num2str(theta_KS_max),'  current solutions: ' , num2str(size(theta_KS_sols.array,2)) ]);
tolerable_KS_direction = Solution_array( IPS_solution, tolerable_CPP_deviations, KS, ...
    theta_KS_max, 'theta_KS_min', 'closest_cp_direction', 'KS', false );
display(['with criterion tolerable_KS_direction = ',num2str(theta_KS_max)] );



%eps_max_solutions = Solution_array( Slip_solution, cpp_deviation_sols, 'eps_ips', eps_max, 'max' ); % for MARESCA&CURTIN RUN
%display(['criterion eps_max=',num2str(eps_max),'  current solutions: ' , num2str(size(eps_max_solutions.array,2)) ]);
eps_max_solutions = Solution_array( IPS_solution, tolerable_CPP_deviations, 'eps_ips', eps_max, 'max' ); 
display(['with criterion eps_max = ',num2str(eps_max)] );



%'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
%theta_NW_max = 10.; % what is allowed, theta_NW_min - what is the smallest angle within the family-transformed_family set
%NW = all_from_family_perms( [1 2 1], false );
% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
%theta_NW_sols = Solution_array( Slip_solution, theta_KS_sols, NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);
%display(['criterion KS_direction_deviation_max=',num2str(theta_KS_max),'  current solutions: ' , num2str(size(theta_NW_sols.array,2)) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering according to Criterion 3: stacking direction
%%% --> with help of - theta_n_max (maximum deviation of block habit-plane to {111}_gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%theta_n_max = 10.; % 2.; % 90.; % % maximum misorientation angle of block habit-plane to {111}_gamma
% specify family near to which habit plane solutions should be searched
% calculation of theta_n - deviation of solution from {111}



% to sort fully reduced solution for most important criterion 
mc_sols = eps_max_solutions.sort( 'theta_h_to_CPP' ) % sort in ascending order for specific property
%det_sols.array(1) % print out best solution

% r = vrrotmat2vec(R)  % r is a four element axis angle rotation vector


