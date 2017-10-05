function [ reduced_sols ] = reduce_Solution_array( initial_sols, prop_string, upper_bound)
% call: [ reduced_sols ] = reduced_solutions( initial_sols, prop_string, prop_interval)
% This function takes an object solution_array and reduces its solutions
% according to the given criteria, upper_bound
% only for the property - 'slip_plane_density' upper_bound actually is a
% lower bound

% in general slip_soluton could also be IPS_solition or something other -
% how to generalize?

switch prop_string
    case 'theta_CPP' % (minimum) misorientation angle of between two close packed planes (cpps) in parent and product phase
        % upper_bound = theta_CPP_max
        reduced_sols = Solution_array( Slip_solution(), initial_sols, cpps_gamma, upper_bound, ...
            'theta_CPP', 'closest_to_cpp', 'cpps_gamma', true);
    case 'theta_h' % (minimum) misorientation angle between habit plane and nearest close packed plane
        % upper_bound = theta_n_max
        reduced_sols = Solution_array( Slip_solution(), initial_sols, ...
                    cpps_gamma, upper_bound, 'theta_h', 'closest_to_h', 'h');
        all_sols.sort( 'theta_n' )
    case 'eps_ips'
        % upper_bound = eps_max
        Solution_array( Slip_solution(), initial_sols, 'eps_ips', upper_bound, 'max' );
        all_sols.sort( 'eps' )
    case 'slip_density'
        % upper_bound = g_min
        reduced_sols = Solution_array( Slip_solution(), initial_sols, 'slip_density', upper_bound, 'min');
    case 'theta_KS_min'
        % upper_bound = thetha_KS_max
        reduced_sols = Solution_array( Slip_solution, initial_sols, KS, theta_KS_max, 'theta_KS_min', 'closest_KS', 'KS', false );
    case 'theta_NW_min'
        % upper_bound = 
        reduced_sols = Solution_array( Slip_solution, tolerable_KS_direction, NW, theta_NW_max, ...
                       'theta_NW_min', 'closest_NW', 'NW', false);
    case 'det'
        % upper_bound = delta_determinant_max
        reduced_sols = Solution_array( Slip_solution, initial_sols, 'det', upper_bound,  det(martensite.U) );
       %%        
        %    case prop_string 
        %        try
        %        catch      
end

% reduced_sols.sort( 'theta_CPP' )

end

