martensite = Martensite(); % creates martensite object
austenite = Base();


austenite.my_base =  [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
martensite.my_base = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
austenite.Bravais_type = 'cubic';
martensite.Bravais_type = 'cubic';
% austenite.Centering = 
% martensite.Centering = 'I';
%%% lattice constants in nm
%%% --> "Modeling the role of external stresses on the 
%%%      austenite-to-bainite phase transformation in 51CrV4 steel",
%%%      M.C. Uslu, D. Canadinc, H.-G. Lambers, S. Tschumak, H.J. Maier,
%%%      Modelling Simul. Mater. Sci. Eng. 19 (2011) 045007, (17pp)
%%%      for 51CrV4: a_gamma = 0.3605   und   a_alpha = 0.2884 in nm
% austenite:  face-centered-cubic (fcc)
a_aust = 3.605 % in Angström
% martensite: body-centered-cubic (bcc)
a_mart = 2.884 % in Angström


austenite.Lp = a_aust*[1 1 1];  
martensite.Lp = a_mart *[1 1 1];

% Correspondance Matrix for B3 Bain
C_am = [0.5,  -0.5,  0.;
        0.5,   0.5,  0.;
        0.,    0.,   1.]; 
    
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

cpps_gamma = all_from_family_perms( [1 1 1] );
% 'misorientation of c.p.p martensite to austenite';
theta_p_max = 2.; % maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice
theta_p_sols = Solution_array( Slip_solution(), all_sols, cpps_gamma, theta_p_max, 'theta_p', 'closest_to_cpp', 'cpps_gamma', true);

% reduce soltuions to ones with g < 20. i.e. at least 20 planes between dislocations
% average number of atom layers before a step due to the (continuum) applied shear occurs (LIS)
g_min = 13.; % could also directly be specified in mod_eigenvalue function e.g. block_symmetric_shear
g_min_sols = Solution_array( Slip_solution(), theta_p_sols, 'g', g_min, 'min'); 

% reduce soltuions to ones with eps < ???. 
eps_max = 0.4;
% Construct reduced array 
eps_max_solutions = Solution_array( Slip_solution(), g_min_sols, 'eps', eps_max, 'max' ); 


theta_n_max = 2.; % maximum misorientation angle of habit-plane to {111}_gamma
% specify family near to which habit plane solutions should be searched
% calculation of theta_n - deviation of solution from {111}
cpp_deviation_sols = Solution_array( Slip_solution(), eps_max_solutions, cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h'); 
% cpp_deviation_sols = Solution_array( Slip_solution(), all_sols, cpps_gamma, theta_n_max, 'theta_a', 'closest_to_a', 'a'); 
% alternatively {557}_gamma could be used here see Iwashita 2011


% specify for which Orientation relations the misorientation should be small in the solutions
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
theta_KS_max = 20.;
KS = all_from_family_perms( [1 1 0], false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )
theta_KS_sols = Solution_array( Slip_solution(), cpp_deviation_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );


% 'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
theta_NW_max = 20.; % what is allowed, theta_NW_min - what is the smallest angle within the family-transformed_family set
NW = all_from_family_perms( [1 2 1], false );
% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
theta_NW_sols = Solution_array( Slip_solution(), theta_KS_sols, NW, theta_NW_max, 'theta_NW_min', 'closest_NW', 'NW', false);

% to sort fully reduced solution for most important criterion 
theta_NW_sols.sort( 'theta_p' ) % sort in ascending order for specific property

% Best solution - determine Habit planes for all symmetry related variants
% for interactions.
% hp_near = martensite.symmetry_related( theta_NW_sols.array(4).h ) 

%% Averaged shape deformation of all "good" solutions (if they are quite similar, which is the case here)

ST_ave = zeros(3);
for i = 1: size( theta_NW_sols.array, 2 )
    ST_ave = ST_ave + abs( theta_NW_sols.array(i).ST );
end
% average all with right determinant...
ST_ave = (ST_ave - theta_NW_sols.array(1).ST - theta_NW_sols.array(4).ST ) / 7

% Round ST_ave and modify the diagonal values of ST_ave such that det(ST)
% is right which yields approximately

        
%delta_detB3_detST_mine = det( ST_mine ) - det(B3);

%ST_s = martensite.symmetry_related( ST_mine );

%for i = 1:size(ST_s,3)
%    strains(:,:,i) = 0.5*(ST_s(:,:,i)' * ST_s(:,:,i) - eye(3) );
%end
%strains