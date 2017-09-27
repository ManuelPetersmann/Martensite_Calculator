function plot_nr_of_solutions( sols, prop_string, upper_bound )
% call: plot_nr_of_solutions_in_property_intervall( sols, prop_string, prop_interval )
% given a solution_array a property string and an intervall/discreteness vector
% [lowerbound,upperbound,point_spacing] for this property makes a plot 

eps_max_solutions = Solution_array( Slip_solution, tolerable_HP_deviations, 'eps_ips', eps_max, 'max' ); 
display(['with criterion eps_max = ',num2str(eps_max)] );

index = 1;
for i = start: delta : stop
    switch prop_string
        case 'theta_p'
            % theta_p_max          
            reduced_sols(index) = count( all_sols, 'theta_p', i);
        case 'theta_n'
            % cpp_max
            reduced_sols(index) = count( all_sols, 'theta_n', i);
        case 'eps'
            % eps_max
            reduced_sols(index) = count( all_sols, 'eps', i);
        case 'g'
            % g_min
            reduced_sols(index) = count( all_sols, 'g', i);
    end
    index = index +1;
end

    function amount = count( all_sols, prop_string, i )
        amount = 0;
        for j = 1 : size(all_sols.array, 2)
            %all_sols.array(j).( prop_string )
            if all_sols.array(j).( prop_string ) < i
                amount = amount+1;
            end 
        end
    end


end

