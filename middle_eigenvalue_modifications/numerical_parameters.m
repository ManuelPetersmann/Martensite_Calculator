%% options for numerical prodecures in "middle_eigenvalue_modifications"

tolerance = 1.e-6; % maximum deviation of lambda_2 to 1
stepwidth_change = 0.5; % smaller and inverse to make increment bigger to escape local minima

%% Khachaturyans slip formulation (g or/= m)   S = 1/g  d \otimes m  % not clear wheter 
% has been replaced by normed shears see below and the g's or m's are
% calculated afterwards
%g_min = 4.;
%g_initial = 100.0; 
%delta_g_initial = 3.; % must be set reasonably so that a solution is found i.e. g_min not to high, delta_g not to high either
% here it reaches as lowest g=5 with delta_g = 5. everything lower will be disregarded

%% slip formulation with unit vectors  S = eps  d_unit \otimes m_unit
eps_max = 1. ; % characteristic shape deformation = shear amplitude in the case of simple shear
eps_initial = 0.;
delta_eps_initial = 0.005;

%% other variables that are used in all solution schemes
I = eye(3);
isol = 0; % counter for number of solutions

