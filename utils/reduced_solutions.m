function [ reduced_sols ] = reduced_solutions( all_sols, start, stop, delta, prop_string)


switch prop_string
    case 'theta_p'
        % theta_p_max
        all_sols.sort( 'theta_p' )
    case 'theta_n'
        % cpp_max
        all_sols.sort( 'theta_n' )
        
    case 'eps'
        % eps_max
        all_sols.sort( 'eps' )
        
    case 'g'
        % g_min
        all_sols.sort( 'g' )
end


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

