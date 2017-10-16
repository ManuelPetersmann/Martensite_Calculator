function [sorted_sols] = sort_solutions( unsrt_sols, prop_name )
% casll sort_solutions(unsrt_sols, prop_name )
% unsrt_sols ...martensite object with sortable properties
% 'eps_ips', 'slip_density', 'theta_h', 'theta_CPP', 'theta_KS_min',
% 'theta_NW_min', 'det'
% this function is only needed if a copy of the sorted sols Solution arrays
% is needed

switch prop_name
    case 'eps_ips'
        sorted_sols = unsrt_sols.sort( 'eps_ips' );
    case 'stepwidth'
        sorted_sols = unsrt_sols.sort( 'stepwidth' );
    case 'theta_h_to_cpp'
        sorted_sols = unsrt_sols.sort( 'theta_h_to_cpp' );
    case 'delta_determinant_max'
        sorted_sols = unsrt_sols.sort( 'delta_determinant_max' );    
    case 'theta_CPPs'
        sorted_sols = unsrt_sols.sort( 'theta_CPPs' );
    case 'theta_KS_min'
        sorted_sols = unsrt_sols.sort( 'theta_KS_min' );
    case 'theta_NW_min'
        sorted_sols = unsrt_sols.sort( 'theta_NW_min' );
    case 'theta_preferred_ILSdir_to_h'
        sorted_sols = unsrt_sols.sort('theta_preferred_ILSdir_to_h');
end
%det_sols.sort( 'dir_of_smallest_def' );
%det_sols.sort( 'dir_of_largest_def' );


end
