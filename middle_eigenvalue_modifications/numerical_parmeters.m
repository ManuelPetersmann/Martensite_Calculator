%% options for numerical prodecures in "middle_eigenvalue_modifications"

solutions = Solution_array( Slip_solution() ); % Construct array with type of solution -> After this line, Solution_array.array is no longer a double 
epsilon = 1.e-12; % accuracy for middle valued eigenvalue
g_min = 4.;
g_initial = 100.0; 
delta_g_initial = 5.; % must be set reasonably so that a solution is found i.e. g_min not to high, delta_g not to high either
% here it reaches as lowest g=5 with delta_g = 5. everything lower will be disregarded
I = eye(3);
isol = 0; % counter for number of solutions

%k = 0.5 % stepwidth scaling factor