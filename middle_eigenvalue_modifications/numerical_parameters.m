%% options for numerical prodecures in "middle_eigenvalue_modifications"

epsilon = 1.e-8; % maximum deviation of lambda_2 to 1
stepwidth_change = 0.5; % smaller and inverse to make increment bigger to escape local minima

%% Khachaturyans slip formulation (g or/= m)   S = 1/g  d \otimes m  % not clear wheter 
g_min = 4.;
g_initial = 100.0; 
delta_g_initial = 1.; % must be set reasonably so that a solution is found i.e. g_min not to high, delta_g not to high either
% here it reaches as lowest g=5 with delta_g = 5. everything lower will be disregarded

%% slip formulation with unit vectors  S = eps  d_unit \otimes m_unit
delta_eps = abs(0.01);


%% other variables that are used in all solution schemes
I = eye(3);
isol = 0; % counter for number of solutions

