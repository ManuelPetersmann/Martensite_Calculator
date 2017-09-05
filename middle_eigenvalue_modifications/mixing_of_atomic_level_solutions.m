function [solutions] = mixing_of_atomic_level_solutions(solutions)
% possible calls: mixing_of_atomic_level_solutions(solutions)
% This function takes deformations F_i of atomic habit plane solutions averages
% them according to   x * F_1 + (1-x) F_2 = F_composite
% and checks again for possible habit plane solutions (lambda_2 = 1)

% Note: THE AVERAGE OF TWO IPS ALWAYS IS AN IPS (except for the small error that is made because the determinant is not invariant to addition)

%% Construct array with type of solution -> After this line, Solution_array.array is no longer a double 
solutions = Solution_array( Slip_solution() ); 

%% set numerical parameters (see file numerical_parameters.m)
numerical_parameters;

% loop over slip system combinations
for is1 = 1:(size(solutions.array,2)-1)
    for is2 = (is1+1):size(solutions.array,2)
        F1 = solutions.array(is1).ST;
        F2 = solutions.array(is2).ST;
       
        dxi = 0.01;
        xi = 0.5; % start at equal fractions of both deformations
        
        F_composite = xi* F1 + (1-xi)*F2;

        %% calculate solution
        % calculate invariant plane vector n_i etc.
        try
            [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F_composite, I, 1.e-3 );
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
        catch
            continue
        end
        
        if 
        
        isol = isol + 2
        
        % Create Slip_solution objects and append them to object array
        solutions.array( isol-1 ) =  Slip_solution(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, eps_s, n, d );
        solutions.array( isol )   =  Slip_solution(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, eps_s, n ,d );
         
    end % end of loop 1
end % end of loop 2


fprintf('number of potential solutions found: n_sol = %i :\n', isol)

end


