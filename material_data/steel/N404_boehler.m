% Zusammensetzung: 15wt% Cr - rest vernachlässigbar...

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
a_aust = 3.595; % 3.6017264; % for 200 Grad Celsius, 
a_mart = 2.885  %2.8807346; % for 200 Grad Celsius, 
austenite.Lp = a_aust*[1 1 1];  
                                              % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
martensite.Lp = a_mart *[1 1 1];  

%% Correspondence Matrix for B3 Bain
% e_mart = C_am * e_aust; infinetively many choices. Take one for slip
% system transformation....

C_am = [0.5,  -0.5,  0.;
        0.5,   0.5,  0.;
        0.,    0.,   1.]; 
%correspondance_matrix(austenite, martensite, m1, m2);

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
% slip systems b.c.c  or  f.c.c
% since the shear is a substantial part of the transformation only 
% shear systems which are favorable in the b.c.c. lattices are considered. 
% the plane and direction families are {110}_alpha, {112}_alpha,
% <111>_alpha, <110>_alpha
plane_families =     [ [1 1 0] ;
                       [1 1 2] ];
direction_families = [ [1 1 1]; 
                       [1 1 0] ];
[ ns, ds ] = independent_slipsystems(  plane_families, direction_families );
% transform them to austenite

% highly symmetric mirror planes from bcc
% {001} family
ms = all_from_family_perms( [0 0 1] );
% {011} family
ms = cat(1, ms, all_from_family_perms( [0 1 1] ) );

%% calculate possible solutions and store solution objects in an object array
all_sols = block_symmetric_doubleshear( B3, cp, ms, ns, ds);
% all_sols = block_symmetric_shear( B3, cp, ms, ns, ds);

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
selection_criteria_maraging;

%% reduce solutions to ones with g < 20. i.e. at least 20 planes between dislocations
% average number of atom layers before a step due to the (continuum) applied shear occurs (LIS)
%g_min = 10.; % 5.; % could also directly be specified in mod_eigenvalue function e.g. block_symmetric_shear
g_min_sols = Solution_array( Slip_solution_doubleshear(), theta_p_sols, 'g', g_min, 'min'); 


%% reduce soltuions to ones with eps < something 
%eps_max = 0.3; % 100.; 
eps_max_solutions = Solution_array( Slip_solution_doubleshear(), g_min_sols, 'eps', eps_max, 'max' ); 


%%
%cpps_gamma = all_from_family_perms( [1 1 1] );
% 'misorientation of c.p.p martensite to austenite';
%theta_p_max = 5.; % 90.;  maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice
theta_p_sols = Solution_array( Slip_solution_doubleshear(), det_sols, cpps_gamma, theta_p_max, 'theta_p', 'closest_to_cpp', 'cpps_gamma', true);


%% 
%theta_n_max = 2.; % 90.; % % maximum misorientation angle of block habit-plane to {111}_gamma
% specify family near to which habit plane solutions should be searched
% calculation of theta_n - deviation of solution from {111}
cpp_deviation_sols = Solution_array( Slip_solution_doubleshear(), eps_max_solutions, cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h'); 
% cpp_deviation_sols = Solution_array( Slip_solution_doubleshear(), all_sols, cpps_gamma, theta_n_max, 'theta_a', 'closest_to_a', 'a'); 
% alternatively {557}_gamma could be used here see Iwashita 2011


%% specify for directions the misorientation should be small in the solutions
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
%theta_KS_max = 3.; % had 6 here in last calculations
%KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )
theta_KS_sols = Solution_array( Slip_solution_doubleshear(), cpp_deviation_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );


% %'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
%theta_NW_max = 10.; % what is allowed, theta_NW_min - what is the smallest angle within the family-transformed_family set
%NW = all_from_family_perms( [1 2 1], false );
% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
theta_NW_sols = Solution_array( Slip_solution_doubleshear(), theta_KS_sols, NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);

%% Added: March 2017
%delta_determinant_max = 0.0001;
det_sols = Solution_array( Slip_solution_doubleshear(), all_sols, 'det', delta_determinant_max,  detB3); 
%
% to sort fully reduced solution for most important criterion 
qi_sols = det_sols.sort( 'theta_p' ); % sort in ascending order for specific property
%theta_NW_sols.array(1) % print out best solution







%% Best solution - determine Habit planes for all symmetry related variants
% for interactions.
% hp must correspond to the solution below for same order!!!
% Vorzeichen meiner Lösung sind so wie die von Lösung Nr4!
%hp_near = martensite.symmetry_related( theta_NW_sols.array(4).h ) 

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
%             0.1100    0.1100    0.8900 ] % 5.4.2017
        
% delta_detB3_detST_mine = det( ST_mine ) - det(B3)

%ST_s = martensite.symmetry_related( ST_mine );

%ST_s = reduce_matrices(ST_s);

%write_strain_from_ST( 'Block_strains_large', ST_s, 1 ); % write strains in largestrain



%% just for check-printing
%for i = 1:size(ST_s,3)
%    strains(:,:,i) = 0.5*(ST_s(:,:,i)' * ST_s(:,:,i) - eye(3) );
%end
%strains


  



       



