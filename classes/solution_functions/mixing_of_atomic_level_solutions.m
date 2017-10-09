function [composite_solutions] = mixing_of_atomic_level_solutions(atomic_solutions, prop_for_tol, tol_prop_doMix) %, minimize_prop )
% call: mixing_of_atomic_level_solutions(atomic_solutions, prop_for_tol, tol_prop_doMix, extremize_prop )
% This function takes deformations F_i of atomic habit plane solutions and averages
% solutions fullfilling compatibility critera given by:
% 'prop_for_tol'... string of property that should be averaged
% tol_prop_doMix... tolerance of property deviation of atomic solutions
% permitting the calculation of an average - following a linear rule of mixture: 
% prop_composite = x * prop_lath_sol1 + (1-x) prop_lath_sol2     e.g.
% F_composite = x*F_lath1 + (1-x)*F_lath2
% with x being the volume fraction of solution 1.
% Note, that for such an average the IPS condition is approximately
% maintained, which has been checked to be true for various combinations.
%
% The variable 'minimize_prop' can be:
%
% 1 ) 'eps': the magnitude of the linearly mixed shape strain vector as a measure of
% self accommodation.  % F_composite = I + eps * d_composite \otimes h_composite
% min_norm2( eps1*d1 + eps2*d2)--> linear constrained optimization
%
% 2 ) | F_composite -I |                   - displacement gradient
% 3 ) | F_composite^T F_composite - I |    - c.f. Green-Lagrange (no rotation)
%
% -)F_composite should be mostly volumetric, i.e.  min|F_composite - % F^spherical| --> nonlinear constrained optimization
% where the diagonal entries of F^spherical are determined from the volume 
% change (det(Bain)) - used e.g. in Qi2014 - though I think this is not reasonable


% Note: THE AVERAGE OF TWO IPS ALWAYS IS AN IPS (except for the small 
% error that is made because the determinant is not invariant to addition)
% In general the linear mixture rule is only valid in the case of
% n->infinity (minor relations) see Bhattacharya - Microstructures of
% martensties - p.131.

%% Construct array with type of solution -> After this line, Solution_array.array is no longer a double 
composite_solutions = Solution_array( Composite_solution() ); 

neglected_due_to_eigenvalue_change_during_mixing = 0;
% loop over slip system combinations
for is1 = 1:(size(atomic_solutions.array,2)-1)
    for is2 = (is1+1):size(atomic_solutions.array,2)
        sol1 = atomic_solutions.array(is1);
        sol2 = atomic_solutions.array(is2);
        
        % do not mix variants not fullfilling predefined criteria
        % e.g. habit plane deviation from {111}_aust or something else
        if abs(sol1.(prop_for_tol) - sol2.(prop_for_tol)) > tol_prop_doMix
            break
        end 
        %if strcmp(minimize_prop,' shape strain vector')
        [x_eps, eps_block] = mixture_vecs_lin_least_squares_opt( cat(2,sol1.eps_ips*sol1.d , sol2.eps_ips*sol2.d) );
        F_comp_eps = reshape( cat(2,reshape(F1,9,1),reshape(F1,9,1)) * x_eps , 3,3) ;
        % it has been found that only for F_comp_eps the determinant is
        % invariant! 
        [x_dis, d_dis] = mixture_matrix_lin_least_squares_opt(cat(3,sol1.ST,sol2.ST));
        % removed third output of above fucntion - ,  F_comp_dis_opt] since
        % determinant is not invariant - same for function call below
        [x_gl,  d_gl] = frob_min_green_lagrange_composite( sol1.ST,sol2.ST, x_eps(1) );
         
        % calculate composite IPS properties of F_comp_eps
        try
            [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F_comp_eps, I, 1.e-3 );
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
        catch
            neglected_due_to_eigenvalue_change_during_mixing = neglected_due_to_eigenvalue_change_during_mixing + 2;
            continue % do not safe solution if
        end
        
        isol = isol + 2;
        % Create Slip_solution objects and append them to object array
        composite_solutions.array( isol-1 ) =  Composite_solution(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, ...
        eps_block, x_eps, F_comp_eps, x_dis, d_dis, x_gl, d_gl );
        composite_solutions.array( isol )   =  Composite_solution(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, ...
        eps_block, x_eps, F_comp_eps, x_dis, d_dis, x_gl, d_gl );
        
    end % end of loop 1
end % end of loop 2

display( ['number of potential solutions found: n_sol = ', num2str(isol) ] )
display( [ num2str(neglected_due_to_eigenvalue_change_during_mixing), ' solutions neglected due to change of eigenvalues'])
end