clear all; close all; clc;

%% Set up coordinate systems and lattice parameters
martensite = Martensite(); % creates martensite object
austenite = Base();

% conventional bases
austenite.my_base =  [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
martensite.my_base = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
austenite.Bravais_type  = 'cubic';
martensite.Bravais_type = 'cubic';
% austenite.Centering = 
martensite.Centering = 'I';
a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 
austenite.Lp = a_aust*[1 1 1];  % 3.5975576 % {1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
                                              % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
martensite.Lp = a_mart *[1 1 1];  % 2.8807346   

%% Correspondence Matrix for B3 Bain
% e_mart = C_am * e_aust; infinetively many choices. Take one for slip
% system transformation....

C_am = [0.5,  -0.5,  0.;
        0.5,   0.5,  0.;
        0.,    0.,   1.]; 
%correspondence_matrix(austenite, martensite, m1, m2);

% define Bain-strain
eta1 = (a_mart/a_aust)*sqrt(2);
eta3 = a_mart / a_aust; % this is one form of three pcossible for the bain strain

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf 1.0 getuned
% werden. Die differenz ist also (n1-1).
B3 = [eta1 0    0   
       0  eta1  0
       0  0  eta3];
martensite.U = B3;

display('Volume change in percent is:');
detB3 = det(B3)

%vars = martensite.variants()

cp = B3*C_am;


%% assemble slip systems in alpha
% slip systems b.c.c  for P^(2)
% following Maresca&Curtin (2017) we consider
% the plane and direction families for P^(2) are {110}_alpha
% and <111>_alpha
plane_families_P2 =     [ [1 1 0] ];
direction_families_P2 = [ [1 1 1] ];
% This has up to now not been considered in any PTMC theory calculation..
[ ns_P2, ds_P2 ] = independent_slipsystems_normal(  plane_families_P2, direction_families_P2 );


% slip systems b.c.c  for P^(3)
% following Maresca&Curtin (2017) we consider
% the plane and direction families for P^(3) are {111}_gamma
% and <101>_gamma
plane_families_P3 =     [ [1 1 1] ];
direction_families_P3 = [ [1 0 1] ];
[ ns_P3, ds_P3 ] = independent_slipsystems_normal(  plane_families_P3, direction_families_P3 );
% % % % transform from fcc to bcc
% % % ns_P3 = (inverse(cp)' * ns_P3')';
% % % ds_P3 = (cp * ds_P3')';

% highly symmetric mirror planes from bcc (not needed for calculations after Maresca&Curtin,
% but we keep them here to test the solutions of our ansatz with use of the slip-systems considered by Maresca&Curtin)
% {001} family
ms = all_from_family_perms( [0 0 1] );
% {011} family
ms = cat(1, ms, all_from_family_perms( [0 1 1] ) ); %This is a more reasonable plane - for a Block boundary!
ms = cat(1, ms, all_from_family_perms( [1 1 2] ) ); % 17.08.2017 - (001)_aust mirror planes transform to (112)_mart 
% for cubic to cubic (for cubic to tetragonal they transform to (110)_mart.


% % ns = cat(1, ns_P2, ns_P3 );
% % ds = cat(1, ds_P2, ds_P3 );

%% calculate possible solutions and store solution objects in an object array
% all_sols = block_symmetric_doubleshear( B3, cp, ms, ns, ds); % for testing of Maresca&Curtins slip systems 
all_sols_specific_slipsys = block_symmetric_doubleshear_specific_slipsys(B3, cp, ms, ns_P2, ds_P2, ns_P3, ds_P3 );
% all_sols = block_symmetric_shear( B3, cp, ms, ns, ds);

%% further checks if solution is appropriate - reduction of total solutions one at a time
% criteria for selection of solutions:
% -) angular deviation to nearest cpp theta_n = min{ angle(n, {111} ) } - must be small
% -) slip density parameter g... number of planes between - steps - must be high
% -) OR's defiations: theta_p = min[ {111}_gamma < {011}_alpha = AL^-T * {111}_gamma } - criterion NR_1 see Qi2013 p.28
%    theta_KS and theta_NW
% -) shape strain - eps_0 - must be small
% -) determinant must be invariant (could change due to additive mixture of matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering according to Criterion 1: plastic accomodation
%%% --> with help of - eps_max (maximum magnitude of shear)
%%% [- g_min (minimal slip plane spacing)  --> not applicable for
%%% MARESCA&CURTIN]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all selection parameters have been moved to a seperate file to enable
% easier comparison between methods - here they have been commented.
% load all parameters from the file:
selection_criteria_maraging;

%% reduce solutions to ones with g < 20. i.e. at least 20 planes between dislocations
% average number of atom layers before a step due to the (continuum) applied shear occurs (LIS)
%g_min = 10.; % 5.; % could also directly be specified in mod_eigenvalue function e.g. block_symmetric_shear
g_min_sols = Solution_array( Slip_solution(), all_sols_specific_slipsys, 'g', g_min, 'min'); 


%% reduce soltuions to ones with eps < something 
%eps_max = 0.3; % 100.; 
eps_max_solutions = Solution_array( Slip_solution(), g_min_sols, 'eps', eps_max, 'max' ); 
% % % % % eps_max_solutions_test = Solution_array( Slip_solution(), all_sols, 'eps', eps_max, 'max' ); 
% % % % % eps_max_solutions_test_crit1 = Solution_array( Slip_solution(), g_min_sols, 'eps', eps_max, 'max' ); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering according to Criterion 2: OR selection
%%% --> with help of - theta_p_max (maximum deviation of gamma- and alpha-cpps)
%%%                  - theta_KS_max (maximum deviation from KS-OR)
%%%              and - theta_NW_max (maximum deviation from NW-OR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cpps_gamma = all_from_family_perms( [1 1 1] );
% 'misorientation of c.p.p martensite to austenite';
%theta_p_max = 5.; % 90.;  maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice
theta_p_sols = Solution_array( Slip_solution(), eps_max_solutions, cpps_gamma, theta_p_max, 'theta_p', 'closest_to_cpp', 'cpps_gamma', true);
% % % % % theta_p_sols_crit2 = Solution_array( Slip_solution(), all_sols, cpps_gamma, theta_p_max, 'theta_p', 'closest_to_cpp', 'cpps_gamma', true);


% specify for directions the misorientation should be small in the solutions
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
%theta_KS_max = 6.; %3.; % had 6 here in last calculations
% KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )
theta_KS_sols = Solution_array( Slip_solution(), theta_p_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );
% % % % % theta_KS_sols_test = Solution_array( Slip_solution(), all_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );
% % % % % theta_KS_sols_crit2 = Solution_array( Slip_solution(), theta_p_sols_crit2, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );

%'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
%theta_NW_max = 10.; % what is allowed, theta_NW_min - what is the smallest angle within the family-transformed_family set
% NW = all_from_family_perms( [1 2 1], false );
% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
theta_NW_sols = Solution_array( Slip_solution(), theta_KS_sols, NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);
% % % % % theta_NW_sols_test = Solution_array( Slip_solution(), all_sols, NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);
% % % % % theta_NW_sols_test_crit2 = Solution_array( Slip_solution(), theta_KS_sols_crit2, NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering according to Criterion 3: stacking direction
%%% --> with help of - theta_n_max (maximum deviation of block habit-plane to {111}_gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%theta_n_max = 10.; % 2.; % 90.; % % maximum misorientation angle of block habit-plane to {111}_gamma
% specify family near to which habit plane solutions should be searched
% calculation of theta_n - deviation of solution from {111}
cpp_deviation_sols = Solution_array( Slip_solution(), theta_NW_sols, cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h'); 
% cpp_deviation_sols = Solution_array( Slip_solution(), all_sols, cpps_gamma, theta_n_max, 'theta_a', 'closest_to_a', 'a'); 
% alternatively {557}_gamma could be used here see Iwashita 2011
% % % % % cpp_deviation_sols_test_crit3 = Solution_array( Slip_solution(), all_sols, cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h'); 


% Added: March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering according to Criterion 4: deviation of det(F)
%%% --> with help of - delta_determinant_max (maximal deviation of det(F) during calculations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%delta_determinant_max = 0.0001;
det_sols = Solution_array( Slip_solution(), cpp_deviation_sols, 'det', delta_determinant_max,  detB3 ); 
% % % % % det_sols_test_crit4 = Solution_array( Slip_solution(), all_sols, 'det', delta_determinant_max,  detB3 ); 
%
% to sort fully reduced solution for most important criterion 
det_sols.sort( 'theta_p' ) % sort in ascending order for specific property
%det_sols.array(1) % print out best solution

% r = vrrotmat2vec(R)  % r is a four element axis angle rotation vector

%% CODE ENDS HERE





% Best solution - determine Habit planes for all symmetry related variants
% for interactions.
% hp must correspond to the solution below for same order!!!
% Vorzeichen meiner Lösung sind so wie die von Lösung Nr4!
% % % hp_near = martensite.symmetry_related( det_sols.array(4).h ) 

%eshelby = martensite.symmetry

%% Averaged shape deformation of all "good" solutions (if they are quite similar, which is the case here)

% ST_ave = zeros(3);
% for i = 1: size( Ni9_theta_NW_sols.array, 2 )
%     ST_ave = ST_ave + abs( Ni9_theta_NW_sols.array(i).ST );
% end
% ST_ave = ST_ave / size( Ni9_theta_NW_sols.array, 2 );

% Round ST_ave and modify the diagonal values of ST_ave such that det(ST)
% is right which yields approximately


% ST_ave = [ 1.0785    0.0784    0.0782
%            0.0911    1.0914    0.0912
%            0.1469    0.1470    0.8534 ]
      
% ST_mine = [ 1.0700    0.0710   -0.0710
%             0.0710    1.0725   -0.0710
%             0.1250    0.1250    0.8800 ]

% ST_mine = [ 1.0650    0.0650   -0.0650
%             0.0650    1.0650   -0.0650
%             0.1100    0.1100    0.8900 ]; % 5.4.2017
        
% delta_detB3_detST_mine = det( ST_mine ) - det(B3)

% ST_s = martensite.symmetry_related( ST_mine );
% 
% ST_s = reduce_matrices(ST_s);

% uncomment if needed
%write_strain_from_ST( 'Block_strains_large', ST_s, 1 ); % write strains in largestrain



%% just for check-printing
%for i = 1:size(ST_s,3)
%    strains(:,:,i) = 0.5*(ST_s(:,:,i)' * ST_s(:,:,i) - eye(3) );
%end
%strains


  



       



