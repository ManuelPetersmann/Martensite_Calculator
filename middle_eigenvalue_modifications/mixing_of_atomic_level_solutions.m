function [solutions] = mixing_of_atomic_level_solutions(atomic_solutions) %, name_of_prop_for_tolerance, tolerance_of_prop_allowing_mix, prop_to_minimize )
% possible calls: mixing_of_atomic_level_solutions(solutions)
% This function takes deformations F_i of atomic habit plane solutions and averages
% solutions fullfilling compatibility critera given by:
% 'name_of_prop_for_tolerance' and 'tolerance_of_prop_allowing_mix'
% according to x * F_1 + (1-x) F_2 = F_composite   with x being the volume
% fraction of solution 1.
% For the averages it is assumed that the IPS condition is approximately
% maintained, which has been checked to be true for various combinations.
% Therefore, an incremental search (varying x) is utilized that optimizes
% the property specified in the variable 'prop_to_minimize': 
% F_composite should be mostly volumetric, i.e.  min|F_composite - F^H| - see e.g. Qi - I think this is not necessary!
% Instead ILS characteristics + 
% the shape strain 'eps' of F_composite = I + eps d \otimes h
% should be as small as possible (without changing the determinant!)

% Note: THE AVERAGE OF TWO IPS ALWAYS IS AN IPS (except for the small 
% error that is made because the determinant is not invariant to addition)

%% Construct array with type of solution -> After this line, Solution_array.array is no longer a double 
composite_solutions = Solution_array( Slip_solution() ); 

%% set numerical parameters (see file numerical_parameters.m)
numerical_parameters;

size(atomic_solutions.array,2)
% loop over slip system combinations
for is1 = 1:(size(atomic_solutions.array,2)-1)
    for is2 = (is1+1):size(atomic_solutions.array,2)
        sol1 = atomic_solutions.array(is1);
        sol2 = atomic_solutions.array(is2);
        F1 = sol1.ST;
        F2 = sol2.ST;
        % shape strains are a measure of self accommodation.
        % per definition of self accommodation: Perfect self 
        eps_ips1 = sol1.eps_ips; 
        eps_ips2 = sol1.eps_ips;  
        
        dxi = 0.01;
        xi = 0.5; % start at equal fractions of both deformations
        
        F_composite = xi* F1 + (1-xi)*F2;
        det(F_composite)
        
  % TODO continue      name_of_prop_for_tolerance, tolerance_of_prop_allowing_mix, prop_to_minimize

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
        
        %if 
        
        isol = isol + 2
        
        % Create Slip_solution objects and append them to object array
        composite_solutions.array( isol-1 ) =  Slip_solution(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, eps_s, n, d );
        composite_solutions.array( isol )   =  Slip_solution(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, eps_s, n ,d );
         
    end % end of loop 1
end % end of loop 2


fprintf('number of potential solutions found: n_sol = %i :\n', isol)

end


