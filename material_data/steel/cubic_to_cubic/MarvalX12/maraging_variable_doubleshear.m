clear all; 
clc;

a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 

Bain_and_Correspondence;

count_directions_extra = true;

% Code ab hier in GUI Code einbauen!


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
direction_families_fcc = [ [1 1 0]; [1 1 2] ];
[austenite.slip_planes, austenite.slip_directions] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);
%[ ns_parent, ds_parent] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);

martensite.considered_plasticity = 3; % 1-mart, 2-aust, 3-both mart and aust slip systems



%% calculate possible solutions and store solution objects in an object array
%martensite.IPS_solutions.array = 
martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite, austenite);



%% further checks if solution is appropriate - reduction of total solutions one at a time
% criteria for selection of solutions:
% -) angular deviation to nearest cpp theta_n = min{ angle(n, {111} ) } - must be small
% -) slip density parameter g... number of planes between - steps - must be high
% -) OR's defiations: theta_CPPs_max = min[ {111}_gamma < {011}_alpha = AL^-T * {111}_gamma } - criterion NR_1 see Qi2013 p.28
%    theta_KS and theta_NW
% -) shape strain - eps_0 - must be small
% -) determinant must be invariant (could change due to additive mixture of matrices

% all selection parameters have been moved to a seperate file to enable
% easier comparison between methods - here they have been commented.
% load all parameters from the file:
selection_criteria_maraging;

% not necessary for laths only for blocks!
det_sols = Solution_array( Slip_solution, martensite.IPS_solutions, 'delta_determinant_max', delta_determinant_max,  det(martensite.U));
display(['with criterion tolerable volume_change_from_averaging = ',num2str(delta_determinant_max)] );




%% ESSENTIAL SELECTION CRITERIA

% Habit plane deviation from experimental observations - theta_h_to_CPP = 20.; %10. % normally between 10 and 20 - see Maresca paper.;  
tolerable_HP_deviations = Solution_array( Slip_solution, det_sols, cpps_gamma, ...
    theta_h_to_CPP, 'theta_h_to_CPP', 'closest_to_h', 'h'); 
display(['with criterion del_habitplane_111gamma_max = ',num2str(theta_h_to_CPP)]);
% alternatively {557}_gamma could be used here see Iwashita 2011


%% 'misorientation of CPP martensite to austenite - planes of OR';
tolerable_CPP_deviations = Solution_array( Slip_solution, tolerable_HP_deviations, cpps_gamma, theta_CPPs_max, ...
    'theta_CPPs', 'closest_to_cpp', 'cpps_gamma', true);
display(['with criterion delta_CPPs_max = ',num2str(theta_CPPs_max)] );
  

%% specify maximum misorientations of solutions to ideal OR directions 
theta_KS_max = 10.; % 6.; 5.; %3.5; % 10.;  had 6 here in last calculations
tolerable_KS_direction = Solution_array( Slip_solution, tolerable_CPP_deviations, KS, ...
    theta_KS_max, 'theta_KS_min', 'closest_cp_direction', 'KS', false );
display(['with criterion tolerable_KS_direction = ',num2str(theta_KS_max)] );





%% EXTENDED selection criteria

%% invariant line crit
% theta_max_ILSdir_to_h = 3.; %
theta_max_ILSdir_to_h = 10.; % 3 reduced it to only 16 from 872 (last reduction step) with criteria from 27.10.17
reduced_sols = Solution_array( Slip_solution, tolerable_KS_direction, austenite.CP_dirs, theta_max_ILSdir_to_h, 'theta_preferred_ILSdir_to_h', 'closest_ILSdir_to_h','KS' ); 
display(['with criterion maximum misorientation of preferred invariant line 110_aust to from invariant habit plane = ',num2str(theta_max_ILSdir_to_h)] );


%% reduce solutions to ones with g < 20. i.e. at least 20 planes between dislocations
% average number of atom layers before a step due to the (continuum) applied shear occurs (LIS)
g_min = 5.; %7.; % 10 % lower bound in numerical_parameters file is 5.!
g_min_sols = Solution_array( Slip_solution, reduced_sols, 'stepwidth', g_min, 'min'); 
display(['with criterion g_min = ',num2str(g_min)] );


%% reduce solutions to ones with eps < something 
eps_max = 0.9;
eps_max_solutions = Solution_array( Slip_solution, g_min_sols, 'eps_ips', eps_max, 'max' ); 
display(['with criterion eps_max = ',num2str(eps_max)] );

%% NW 
% tolerable_NW_direction = Solution_array( Slip_solution, tolerable_KS_direction, NW, theta_NW_max, ...
%     'theta_NW_min', 'closest_NW', 'NW', false);
% display(['with criterion tolerable delta_CPP_max = ',num2str(theta_CPPs_max)] );





%% plots for IPS solution sensitivity - start after reduction to ESSENTIAL SELECTION CRITERIA

% nr = 0;
% for i = 1:length( tolerable_HP_deviations.array)
%     if tolerable_HP_deviations.array(i).axis_angle_rotvec_inclusion(4) < hp_tol
%         nr = nr +1;
%     end
% end
% nr

%
plot_nr_of_solutions_in_property_intervall( reduced_sols, 'eps_ips', [0, 1.0, 41] )

%
% alternatively {557}_gamma could be used here see Iwashita 2011
plot_nr_of_solutions_in_property_intervall( reduced_sols, 'theta_h_to_CPP', [0, 20, 41] )

%
plot_nr_of_solutions_in_property_intervall( reduced_sols, 'theta_preferred_ILSdir_to_h', [0, 10, 41] )

%
plot_nr_of_solutions_in_property_intervall( reduced_sols, 'stepwidth', [0, 20, 41] )

%
plot_nr_of_solutions_in_property_intervall( reduced_sols, 'theta_KS_min', [0, 3, 41] )



%% Post Processing of Lath Solutions

%% sorting
% to sort fully reduced solution for most important criterion 
% mar_sols = det_sols.sort( 'theta_CPP' ); % sort in ascending order for specific property
% theta_NW_sols.array(1) % print out best solution
%sorted_sols = reduced_sols.sort( 'stepwidth' );   


%% slip ambiguity reduction
d1_tol = 1.e-2;
gelockerte_lath_constraints = multiplicity_check_due_to_slip( reduced_sols ); % angle_tol, def_tol  %, d1_tol ); 

% copy last state of solutions except for reduced ones in last step:
glc = reduced_sols; % tolerable_KS_direction;
glc.array = gelockerte_lath_constraints;



%% BLOCK calculations
block_solutions = Solution_array_composite();

% other block contraints...
% theta_hps = 5; % does not reduce anything after so many restrictions were placed on laths...
% theta_intersec_cpdir = 7. % 7 - leads to 35 sols % 6. - leads to zero solutions... 
% block_solutions.mixing_tolerances('theta_intersec_cpdir') = theta_intersec_cpdir;
% block_solutions.mixing_tolerances('theta_hps') = theta_hps;

% lambda2_tol = 1.e-3; % if taken as 1.e-5 than no solutions get sorted out that way
% cof_tol = 1.e-6;
% det_tol = 1.e-6;

% function to investigate solution space of blocks
block_tests(glc, block_solutions, martensite.U); %, lambda2_tol, cof_tol, det_tol);
% starting from: tolerable_KS_direction, i.e. only CP, KS and CPP constraints -> leads to 347 entries - 48 solutions in the
% lambda2_tol = cof_tol = lambda2_tol = 5.e-3
% + crit ILS dir = 6 -> 308 entries,  40 solutions in the end 
% + crit eps_max = 0.9 -> 260 entries, 24 solutions in the end
% + crit g_min > 10  --> 64 entries, 0 solutions

block_solutions = mixing_of_atomic_level_solutions(glc, block_solutions, martensite.U);



  



       



