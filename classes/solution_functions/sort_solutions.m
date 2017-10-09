function [sorted_sols] = sort_solutions( unsrt_sols, prop_name )
% casll sort_solutions(unsrt_sols, prop_name )
% unsrt_sols ...martensite object with sortable properties
% 'eps_ips', 'slip_density', 'theta_h', 'theta_CPP', 'theta_KS_min',
% 'theta_NW_min', 'det'

switch prop_name
    case 'eps_ips'
        sorted_sols = unsrt_sols.sort( 'eps_ips' );
    case 'slip_density'
        sorted_sols = unsrt_sols.sort( 'slip_density' );
    case 'theta_h'
        sorted_sols = unsrt_sols.sort( 'theta_h' );
    case 'det'
        sorted_sols = unsrt_sols.sort( 'det' );    
    case 'theta_CPP'
        sorted_sols = unsrt_sols.sort( 'theta_CPP' );
    case 'theta_KS_min'
        sorted_sols = unsrt_sols.sort( 'theta_KS_min' );
    case 'theta_NW_min'
        sorted_sols = unsrt_sols.sort( 'theta_NW_min' );
end
%det_sols.sort( 'dir_of_smallest_def' );
%det_sols.sort( 'dir_of_largest_def' );


end

