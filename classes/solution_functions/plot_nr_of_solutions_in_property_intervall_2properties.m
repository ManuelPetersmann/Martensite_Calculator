function plot_nr_of_solutions_in_property_intervall_2properties( all_sols, prop_strings, intervalls) %, xlab, ylab )
% call: plot_nr_of_solutions_in_property_intervall( sols, prop_string, intervall )
% given a solution_array a property string and an intervall/discreteness vector
% intervall = [ lowerbound, upperbound, data_points ] cumulative plot 

% eps_max_solutions = Solution_array( Slip_solution, tolerable_HP_deviations, 'eps_ips', eps_max, 'max' ); 
% display(['with criterion eps_max = ',num2str(eps_max)] );

x = linspace( intervalls(1,1), intervalls(1,2), intervalls(1,3) );
y = linspace( intervalls(2,1), intervalls(2,2), intervalls(2,3) );
cumulative = [];
for i = 1:length(x)
    cumulative_x = get_x( all_sols, prop_strings{1}, x(i) );
    if ~isempty( cumulative_x(1).F1 )
        for ii = 1:length(y)
            %length(cumulative_x)
            cumulative(i,ii) = count( cumulative_x, prop_strings{2}, y(ii) );
        end
    end
end


%figure;
%surf(x,y,cumulative) % this was wrong!

figure;
surf(y,x,cumulative)
% xlabel( prop_string,'Interpreter', 'none' ) % x-axis label
% ylabel( 'Nr of IPS solutions' ) % y-axis label
% ax = gca;
%ax.YLim = [0,700];

% NameArray = {'Marker','Tag'};
% ValueArray = {'o','Decaying Exponential';...
%    'square','Growing Exponential';...
%    '*','Steady State'};
% set(S,NameArray,ValueArray)
    
% if nargin > 3
%     xlabel( xlab ) % x-axis label
%     ylabel( ylab ) % y-axis label
% end

   function xx = get_x( all_sols, prop_string, maxi )
        amount = 1;
        xx = IPS_solution;
        if isKey(all_sols.array(1).added_props,  prop_string )
            for j = 1 : size(all_sols.array, 2)
                if all_sols.array(j).added_props( prop_string ) < maxi
                    xx(amount) = all_sols.array(j);
                    amount = amount+1;
                end
            end
        else
            if isprop( all_sols.array(1), prop_string) || isprop( all_sols.array(1).slip, prop_string)
                for j = 1 : size(all_sols.array, 2)
                    if strcmp(prop_string,'max_eps_s') % two values - sort for the smaller one 1/stepwidth \propto eps_s (plastic shear magnitude)
                        if  min(all_sols.array(j).slip.(prop_string)) < maxi
                            xx(amount) = all_sols.array(j);
                            amount = amount+1;
                        end
                    else
                        if all_sols.array(j).( prop_string ) < maxi
                            xx(amount) = all_sols.array(j);
                            amount = amount+1;
                        end
                    end
                end
            else
                error('prop_string is not the name of any property of this object')
            end
        end
    end
%%

    function amount = count( all_sols, prop_string, maxi )
        amount = 0;
        if isKey(all_sols(1).added_props,  prop_string )
            for j = 1 : size(all_sols, 2)
                if all_sols(j).added_props( prop_string ) < maxi
                    amount = amount+1;
                end
            end
        else
            if isprop( all_sols(1), prop_string) || isprop( all_sols(1).slip, prop_string)
                for j = 1 : size(all_sols, 2)
                    if strcmp(prop_string,'max_eps_s') % two values - sort for the smaller one 1/stepwidth \propto eps_s (plastic shear magnitude)
                        if  min(all_sols(j).slip.(prop_string)) < maxi
                            amount = amount+1;
                        end
                    else
                        if all_sols(j).( prop_string ) < maxi
                            amount = amount+1;
                        end
                    end
                end
            else
                error('prop_string is not the name of any property of this object')
            end
        end
    end


end



%         switch prop_string
%         case 'theta_p'
%             % theta_p_max          
%             reduced_sols(index) = count( all_sols, 'theta_p', i);
%         case 'theta_n'
%             % cpp_max
%             reduced_sols(index) = count( all_sols, 'theta_n', i);
%         case 'eps'
%             % eps_max
%             reduced_sols(index) = count( all_sols, 'eps', i);
%         case 'g'
%             % g_min
%             reduced_sols(index) = count( all_sols, 'g', i);