clear all; 
clc;
%format long

a_aust = 3.62; % Consider Temperature dependence? -> see joint file!
a_mart = 2.872; % Consider Temperature dependence? -> see joint file!

Bain_and_Correspondence;

% assemble slip systems in alpha
% since the shear is a substantial part of the transformation only 
% shear systems which are favorable in the b.c.c. lattices are considered. 
% the plane and direction families are {110}_alpha, {112}_alpha,
% <111>_alpha, <110>_alpha
plane_families_bcc =     [ [1 1 0]
                           [1 1 2] ];   % must be written with linebreak or ";" between vectors!                     
direction_families_bcc = [ [1 1 1]
                           [1 1 0] ];
                       
count_directions_extra = true;
                       
% find all possible combination (including different shear directions)
[martensite.slip_planes, martensite.slip_directions] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );
%[ ns_product, ds_product ] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );

plane_families_fcc =     [ [1 1 1] ; [1 1 2] ]; % PET 09.03.18 - added plane 112 - Kelly 1992
direction_families_fcc = [ [1 1 0] ; [1 1 2] ]; % 112 is twinning dislocation... - 112 is a partial! do not take it!
[austenite.slip_planes, austenite.slip_directions] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);
%[ ns_parent, ds_parent] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);

martensite.considered_plasticity = 3; % 1-mart, 2-aust, 3-both mart and aust slip systems


% further checks if solution is appropriate - reduction of total solutions one at a time
% criteria for selection of solutions:
% -) angular deviation to nearest cpp theta_n = min{ angle(n, {111} ) } - must be small
% -) slip density parameter g... number of planes between - steps - must be high
% -) OR's defiations: theta_CPPs_max = min[ {111}_gamma < {011}_alpha = AL^-T * {111}_gamma } - criterion NR_1 see Qi2013 p.28
%    theta_KS and theta_NW
% -) shape strain - eps_0 - must be small
% -) determinant must be invariant (could change due to additive mixture of matrices
%%
% all selection parameters have been moved to a seperate file to enable
% easier comparison between methods - here they have been commented.
% load all parameters from the file:
selection_criteria_mcl_bainite;


%% IPS lath solutions
% specify wheter after each incremental slip the crystal lattice is also
% rotated and the slip systems (large strain crystal plasticity)
update_correspondence = false;

%% calculate possible IPS solutions and store solution objects in an object array
%martensite.IPS_solutions.array = 
martensite.IPS_solutions = doubleshear_variable_shear_mags( martensite, austenite, update_correspondence);


%% not necessary for laths only for blocks!
det_sols = Solution_array( IPS_solution, martensite.IPS_solutions, 'delta_determinant_max', delta_determinant_max,  det(martensite.U));
display(['with criterion tolerable volume_change_from_averaging = ',num2str(delta_determinant_max)] );


%% ESSENTIAL SELECTION CRITERIA

% Habit plane deviation from experimental observations - theta_h_to_CPP = 20.; %10. % normally between 10 and 20 - see Maresca paper.;  
tolerable_HP_deviations = Solution_array( IPS_solution, det_sols, cpps_gamma, ...
    theta_h_to_CPP, 'theta_h_to_CPP', 'closest_to_h', 'h'); 
display(['with criterion del_habitplane_111gamma_max = ',num2str(theta_h_to_CPP)]);
% alternatively {557}_gamma could be used here see Iwashita 2011


%% 'misorientation of CPP martensite to austenite - planes of OR';
tolerable_CPP_deviations = Solution_array( IPS_solution, tolerable_HP_deviations, cpps_gamma, theta_CPPs_max, ...
    'theta_CPPs', 'closest_to_cpp', 'cpps_gamma', true);
display(['with criterion delta_CPPs_max = ',num2str(theta_CPPs_max)] );
  

%% specify maximum misorientations of solutions to ideal OR directions 
%theta_KS_max = 5.26; % 6.; 5.; %3.5; % 10.;  had 6 here in last calculations
tolerable_KS_direction = Solution_array( IPS_solution, tolerable_CPP_deviations, KS, ...
    theta_KS_max, 'theta_KS_min', 'closest_cp_direction', 'KS', false );
display(['with criterion tolerable_KS_direction = ',num2str(theta_KS_max)] );



%% EXTENDED selection criteria

%% invariant line crit
% theta_max_ILSdir_to_h = 3.; %
%theta_max_ILSdir_to_h = 10.; % 3 reduced it to only 16 from 872 (last reduction step) with criteria from 27.10.17
max_KS_dirs_sols = Solution_array( IPS_solution, tolerable_KS_direction, austenite.CP_dirs, theta_max_ILSdir_to_h, 'theta_preferred_ILSdir_to_h', 'closest_ILSdir_to_h','KS' ); 
display(['with criterion maximum misorientation of preferred invariant line 110_aust to from invariant habit plane = ',num2str(theta_max_ILSdir_to_h)] );


%% reduce solutions to ones with g < 20. i.e. at least 20 planes between dislocations
% average number of atom layers between straight dislocation lines in the (continuum) applied shear 
%g_min = 5.; %7.; % 10 % lower bound in numerical_parameters file is 5.!
max_eps_s_sols = Solution_array( IPS_solution, max_KS_dirs_sols, 'max_eps_s', eps_s_max, 'max_shear');   % PET 01.03.18: exchanged 'stepwidth' with 'max_eps_s'
display(['with criterion max_eps_s = ',num2str(eps_s_max)] );


max_eps_s_sum_sols = Solution_array( IPS_solution, max_eps_s_sols, 'eps_s_sum', eps_s_sum_max, 'shear_sum');   % PET 01.03.18: exchanged 'stepwidth' with 'max_eps_s'
display(['with criterion max_eps_s = ',num2str(eps_s_max)] );


%% reduce solutions to ones with eps < something 
%eps_IPS_max = 0.9;
eps_IPS_max_solutions = Solution_array( IPS_solution, max_eps_s_sum_sols, 'eps_ips', eps_ips_max, 'max' ); 
display(['with criterion eps_ips_max = ',num2str(eps_ips_max)] );

%% reduce lath solutions to ones with rotation angle < 5 degree
% red_sols_ILS =    Solution_array( IPS_solution, eps_IPS_max_solutions, 'rotangle_inclusion', max_rot_angle_lath);
% disp([' for a maximum lath rotation angle  < ',num2str(max_rot_angle_lath),'Â°'] );

%% NW 
% tolerable_NW_direction = Solution_array( IPS_solution, tolerable_KS_direction, NW, theta_NW_max, ...
%     'theta_NW_min', 'closest_NW', 'NW', false);
% display(['with criterion tolerable delta_CPP_max = ',num2str(theta_CPPs_max)] );

reduced_sols = eps_IPS_max_solutions;

%% plots for IPS solution sensitivity - start after reduction to ESSENTIAL SELECTION CRITERIA

% nr = 0;
% for i = 1:length( tolerable_HP_deviations.array)
%     if tolerable_HP_deviations.array(i).axis_angle_rotvec_inclusion(4) < hp_tol
%         nr = nr +1;
%     end
% end
% nr

% plot_nr_of_solutions_in_property_intervall( reduced_sols, 'eps_ips', [0, eps_ips_max, 100] )
% 
% % alternatively {557}_gamma could be used here see Iwashita 2011
%plot_nr_of_solutions_in_property_intervall( reduced_sols, 'theta_h_to_CPP', [0, theta_h_to_CPP, 100] )
% 
%plot_nr_of_solutions_in_property_intervall( reduced_sols, 'theta_preferred_ILSdir_to_h', [0, theta_max_ILSdir_to_h, 100] )
% 
% plot_nr_of_solutions_in_property_intervall( reduced_sols, 'max_eps_s', [0.1, eps_s_max, 100] )
% 
% plot_nr_of_solutions_in_property_intervall( reduced_sols, 'eps_s_sum', [0.1, eps_s_sum_max, 100] )
% 
% plot_nr_of_solutions_in_property_intervall( reduced_sols, 'theta_KS_min', [0, theta_KS_max, 100] )
% 
%plot_nr_of_solutions_in_property_intervall( reduced_sols, 'theta_CPPs', [0, theta_CPPs_max, 100] )
% 
% %%
% plot_nr_of_solutions_in_property_intervall( reduced_sols, 'rotangle_inclusion', [-20, 20, 100] ) 
% hold on 
% plot_nr_of_solutions_in_property_intervall( reduced_sols, 'rotangle_inclusion', [0, 20, 100] ) 
% 
% %%
plot_nr_of_solutions_in_property_intervall_2properties( reduced_sols, {'eps_ips','theta_CPPs'}, cat(1,[0, eps_ips_max, 41],[0, theta_CPPs_max, 41]) )
% %%
plot_nr_of_solutions_in_property_intervall_2properties( reduced_sols, {'rotangle_inclusion','eps_ips'}, cat(1,[0, 15, 41],[0, eps_ips_max, 41]) )


%% Post Processing of Lath Solutions

%% sorting
% to sort fully reduced solution for most important criterion 
%mar_sols = reduced_sols.sort( 'eps_ips' ); % sort in ascending order for specific property
% theta_NW_sols.array(1) % print out best solution
%sorted_sols = reduced_sols.sort( 'stepwidth' );   


%% slip ambiguity reduction
% d1_tol = 1.e-2;
reduced_sols = multiplicity_check_due_to_slip( reduced_sols ); % angle_tol, def_tol  %, d1_tol ); 
% 
% % copy last state of solutions except for reduced ones in last step:
% glc = reduced_sols; % tolerable_KS_direction;
% glc.array = gelockerte_lath_constraints;


%% BLOCK calculations... 
% Lass ma das Mal... hab alles rausgelöscht bis auf den letzten
% Kommentierten Code. Weis noch dass ich da mit dem Richard Jurisitz diskutiert habe
% Wer lust hat kann "added_mass" googleln und sich fragen wie man auf sowas
% kommt XD...

 
% [~ , R_total] = polardecomposition( F_tot );
% [angle_inclusion, ax_inclusion] = rotmat_to_axis_angle( R_total );
% added_mass_angle = angle_ax_u + angle_inclusion; % added_mass_angle XD
% if (angle_inclusion < max_rot_angle_inclusion)
%     %if (added_mass_angle < added_mass_angle_tolerance)
%     %[angle_lattice, ax_lattice] = rotmat_to_axis_angle( R_mod );
%     
%     [bool, lambda2] = is_rank_one_connected(F_tot,eye(3), lambda2_ips_tolerance_lath);
%     if bool
%         
%         LT = R_mod*martensite.U;
%         
%         [ theta_CP, closest_111aust_to_CP ] = min_misorientation( austenite.CPPs, LT, true );
%         %theta_CP
%         if theta_CP < theta_CP_max


  



       



