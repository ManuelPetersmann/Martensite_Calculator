function [ reduced_sols ] = reduce_Solution_array( initial_sols, austenite, prop_string, upper_bound)

% TODO TODO TODO not used for now...

% call: [ reduced_sols ] = reduced_solutions( initial_sols, prop_string, prop_interval)
% This function takes an object solution_array and reduces its solutions
% according to the given criteria, upper_bound
% only for the property - 'slip_plane_density' upper_bound actually is a
% lower bound

% in general slip_soluton could also be IPS_solition or something other -
% how to generalize?
metaclass = ?initial_sols;


switch prop_string
    case 'theta_CPPs' % (minimum) misorientation angle of between two close packed planes (cpps) in parent and product phase
        % upper_bound = theta_CPP_max
        reduced_sols = Solution_array( Slip_solution(), initial_sols, austenite.CPPs, upper_bound, ...
            'theta_CPPs', 'closest_CPPs', 'cpps_gamma', true);
    %    Log_sol_info(red_sols,[' for misorietation of CPPs  < ',num2str(theta_CPPs_max),'°'] );
    case 'theta_h_to_CPP' % (minimum) misorientation angle between habit plane and nearest close packed plane
        % upper_bound = theta_h_to_CPP
        reduced_sols = Solution_array( Slip_solution(), initial_sols, ...
            austenite.CPPs, upper_bound, 'theta_h_to_CPP', 'closest_h_to_CPP', 'h');
    %    Log_sol_info(red_sols,[' for habit plane misorientation to CP-planes  < ',num2str(theta_h_to_cpp),'°'] ); 
    case 'eps_ips'
        % upper_bound = eps_max
        Solution_array( Slip_solution(), initial_sols, 'eps_ips', upper_bound, 'max' );
    %    Log_sol_info(red_sols,[' for a shape strain  < ',num2str(eps_ips_max)] );
    case 'stepwidth'
        % upper_bound = g_min
        reduced_sols = Solution_array( Slip_solution(), initial_sols, 'stepwidth', upper_bound, 'min');
    %    Log_sol_info(red_sols,[' for a stepwidth > ',num2str(stepwidth)] );
    case 'theta_KS_min'
        % upper_bound = thetha_KS_max
        reduced_sols = Solution_array( Slip_solution, initial_sols, austenite.CP_dirs, upper_bound, 'theta_KS_min', 'closest_KS', 'KS', false );
    %    Log_sol_info(red_sols,[' for a maximum deviation angle of KS-directions  < ',num2str(theta_KS_max),'°'] );
    case 'theta_NW_min'
        % upper_bound = theta_NW_max
        reduced_sols = Solution_array( Slip_solution, tolerable_KS_direction, NW, upper_bound, ...
            'theta_NW_min', 'closest_NW', 'NW', false);
    %    Log_sol_info(red_sols,[' for a maximum deviation angle of NW-directions  < ',num2str(theta_NW_max),'°'] );
    case 'delta_determinant_max'
        % upper_bound = delta_determinant_max
        reduced_sols = Solution_array( Slip_solution, initial_sols, 'delta_determinant_max', upper_bound,  det(martensite.U) );
    %    Log_sol_info(red_sols,[' for (non-physical) volume change  < ',num2str(delta_determinant_max)] );
    case 'theta_max_ILSdir_to_h'
        % upper_bound = theta_e1_cpp_dir
        reduced_sols = Solution_array( Slip_solution, initial_sols, austenite.CP_dirs, upper_bound, 'theta_preferred_ILSdir_to_h', 'closest_ILSdir_to_h' );
    %    Log_sol_info(red_sols,[' for a maximum deviation angle of preferred invariant line from invariant habit plane < ',num2str(theta_max_ILSdir_to_h)] );
        % new with 6 arguments! - no property is added dynamically this way
        % wrong like this: reduced_sols = Solution_array( Slip_solution, initial_sols, austenite.CP_dirs, upper_bound, 'theta_e1_cp_dir', 'closest_cp_dir_to_e1' );
end

    function Log_sol_info( red_sols,crit )
        if (size( red_sols.array, 2)==1) && isempty(red_sols.array(1).F1)
            disp('No Solution fullfilling specified criteria');
        else
            disp(['Solutions reduced to : ', num2str(length(red_sols.array)), crit ] );
        end
    end

end

