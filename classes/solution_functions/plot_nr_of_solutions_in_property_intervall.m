function plot_nr_of_solutions_in_property_intervall( all_sols, prop_string, intervall) %, xlab, ylab )
% call: plot_nr_of_solutions_in_property_intervall( sols, prop_string, intervall )
% given a solution_array a property string and an intervall/discreteness vector
% intervall = [ lowerbound, upperbound, data_points ] cumulative plot 

% eps_max_solutions = Solution_array( Slip_solution, tolerable_HP_deviations, 'eps_ips', eps_max, 'max' ); 
% display(['with criterion eps_max = ',num2str(eps_max)] );


index = 1;
x = linspace( intervall(1), intervall(2), intervall(3) );
for i = 1:length(x)
    cumulative(index) = count( all_sols, prop_string, x(i) );
    index = index +1;
end

figure;
plot(x,cumulative,'r-x')
xlabel( prop_string,'Interpreter', 'tex' ) % x-axis label
ylabel( 'Nr of IPS solutions' ) % y-axis label
ax = gca; % get-current-axis
%ax.YLim = [0,700];

%NameArray = {'Marker','LineWidth'}; %,'Tag'};
%ValueArray = {'x',2}; % 'Decaying Exponential';...
%    'square','Growing Exponential';...
%    '*','Steady State'};
%set(S,NameArray,ValueArray)

% if nargin > 3
%     xlabel( xlab ) % x-axis label
%     ylabel( ylab ) % y-axis label
% end



    function amount = count( all_sols, prop_string, maxi )
        amount = 0;
        if isKey(all_sols.array(1).added_props,  prop_string )
            for j = 1 : size(all_sols.array, 2)
                if all_sols.array(j).added_props( prop_string ) < maxi
                    amount = amount+1;
                end
            end
        else
            if isprop( all_sols.array(1), prop_string) || isprop( all_sols.array(1).slip, prop_string)
                for j = 1 : size(all_sols.array, 2)
                    switch prop_string
                        case 'max_eps_s' % if strcmp(prop_string,'max_eps_s') % two values - sort for the smaller one 1/stepwidth \propto eps_s (plastic shear magnitude)
                            if  all_sols.array(j).slip.(prop_string) < maxi
                                amount = amount+1;
                            end
                        case 'eps_s_sum' %else
                            if  all_sols.array(j).slip.(prop_string) < maxi
                                amount = amount+1;
                            end
                        otherwise
                            if all_sols.array(j).( prop_string ) < maxi
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