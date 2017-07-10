martensite = Martensite(); % creates martensite object
austenite = Base();

austenite.my_base =  [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
martensite.my_base = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.];
austenite.Bravais_type = 'cubic';
martensite.Bravais_type = 'cubic';
% austenite.Centering = 
% martensite.Centering = 'I';
a_aust = 3.56; % see Qi et al. 2014
a_mart = 2.869; % see Qi et al. 2014
austenite.Lp = a_aust*[1 1 1];  % 3.5975576 %{1.0, 1.0, 1.0,  pi/2,   pi/2,   pi/2; ...
           % 'a'  'b'  'c'  'alpha'  'beta'  'gamma'}
martensite.Lp = a_mart *[1 1 1];  % 2.8807346     

%% Correspondance Matrix for B3 Bain
C_am = [0.5,  -0.5,  0.;
        0.5,   0.5,  0.;
        0.,    0.,   1.]; 
% % Correspondance Matrix for B3 Bain (other basis)
% C_am = [ 0.5,   0.5,  0.;
%         -0.5,   0.5,  0.;
%          0.,    0.,   1.];       
%correspondance_matrix(austenite, martensite, m1, m2);

% define Bain-strain
eta1 = (a_mart/a_aust)*sqrt(2);
eta3 = a_mart / a_aust; % this is one form of three possible for the bain strain

% temprary: rounding of values for direct comparison with results from paper
% n1 = round ( ((a_mart/a_aust)*sqrt(2))*1000 )/1000;
% n3 = round ( (a_mart / a_aust)*1000 )/1000;

% Der mittlere Eigenwert ist hier also n1. Dieser soll auf eins getuned
% werden. Die differenz ist also (n1-1).
B3 = [eta1 0    0   
       0  eta1  0
       0  0  eta3];
martensite.U = B3;

det_B3 = det(B3) % volume change
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

% to reproduce solution 1 of solutions in paper
% mss = [[0 1 0]];
% nss = [[1 1 0]; [0 1 -1]];
% dss = [[1 -1 0]; [0 -1 -1]];
% to reproduce solution 2 of solutions in paper
% mss = [[1 0 1]];
% nss = [[1 -1 0]; [0 1 -1]];
% dss = [[1 1 0]; [0 -1 -1]];

%% calculate possible solutions and store solution objects in an object array
all_sols = block_symmetric_doubleshear( B3, cp, ms, ns, ds);
%all_sols = block_symmetric_shear( B3, cp, ms, ns, ds);

%% further checks if solution is appropriate - reduction of total solutions one at a time
% criteria for selection of solutions:
% -) angular deviation to nearest cpp theta_n = min{ angle(n, {111} ) } - must be small
% -) slip density parameter g... number of planes between - steps - must be high
% -) OR's defiations: theta_p = min[ {111}_gamma < {011}_alpha = AL^-T *{111}_gamma } - criterion NR_1 see Qi2013 p.28
%    theta_KS and theta_NW
% -) shape strain - eps_0 - must be small


cpps_gamma = all_from_family_perms( [1 1 1] );
% 'misorientation of c.p.p martensite to austenite';
theta_p_max = 2. % maximum misorientation angle of cpps gamma & alpha - due to Qi,Khachaturyan 2013
% misorientation-angle theta_p between the closed-packed planes (cpp) of alpha {110} and gamma {111} lattice
Ni9_theta_p_sols = Solution_array( Slip_solution(), all_sols, cpps_gamma, theta_p_max, 'theta_p', 'closest_to_cpp', 'cpps_gamma', true);

% reduce soltuions to ones with g < 20. i.e. at least 20 planes between dislocations
% average number of atom layers before a step due to the (continuum) applied shear occurs (LIS)
g_min = 13.; 
Ni9_g_max_sols = Solution_array( Slip_solution(), Ni9_theta_p_sols, 'g', g_min, 'min' ); 

% reduce soltuions to ones with eps < eps_max. 
eps_max = 0.4;
% Construct reduced array 
eps_max_solutions = Solution_array( Slip_solution(), Ni9_g_max_sols, 'eps', eps_max, 'max' ); 

theta_n_max = 2. % maximum misorientation angle of habit-plane to {111}_gamma
% specify family near to which habit plane solutions should be searched
% calculation of theta_n - deviation of solution from {111}
Ni9_cpp_deviation_sols = Solution_array( Slip_solution(), eps_max_solutions, cpps_gamma, theta_n_max, 'theta_n', 'closest_to_h', 'h'); 
% cpp_deviation_sols = Solution_array( Slip_solution(), all_sols, cpps_gamma, theta_n_max, 'theta_a', 'closest_to_a', 'a'); 
% alternatively {557}_gamma could be used here see Iwashita 2011

% specify for which Orientation relations the misorientation should be small in the solutions
% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
theta_KS_max = 90.;
KS = all_from_family_perms( [1 1 0] ); %, false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
% omega in Paper = theta_KS here; angle( [1 -1 0]_gamma, [1-10]_alpha )
Ni9_theta_KS_sols = Solution_array( Slip_solution(), Ni9_cpp_deviation_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );


% 'Nishiyama Wassermann directions: [112]_aust || [110]_mart or equivalently [112]_aust || [110]_mart';
theta_NW_max = 90.; % 
NW = all_from_family_perms( [1 2 1] ); %, false );
% omega - 5.26 in Paper = theta_NW; angle( [1 -2 1], [1 0 -1] )
Ni9_theta_NW_sols = Solution_array( Slip_solution(), Ni9_theta_KS_sols, NW, theta_NW_max, 'theta_NW_min', 'closest_to_NW', 'NW', false);

% to sort fully reduced solution for most important criterion 
Ni9_theta_NW_sols.sort( 'theta_p' ); % sort in ascending order for specific property

% 'looooooooooooooooooooooooooooooooooooooooooooooooooooooooooool'
% [ theta_NW_min, closest_to_NW ] = min_misorientation( NW, Ni9_theta_NW_sols.array(1).LT, false )
%  
% [ theta_KS_min, closest_to_KS ] = min_misorientation( KS, Ni9_theta_NW_sols.array(1).LT, false )

%% Averaged shape deformation of all "good" solutions (if they are quite similar, which is the case here)

% ST_ave = zeros(3);
% for i = 1: size( Ni9_theta_NW_sols.array, 2 )
%     ST_ave = ST_ave + abs( Ni9_theta_NW_sols.array(i).ST );
% end
% ST_ave = ST_ave / size( Ni9_theta_NW_sols.array, 2 );

% Round ST_ave and modify the diagonal values of ST_ave such that det(ST)
% is right which yields approximately - ST_mine

ST_min = [1.1000    0.1000   -0.1000
          0.1000    1.1000   -0.1000
          0.1500    0.1500    0.8455];
    
vars = martensite.symmetry_variants( ST_ave )






