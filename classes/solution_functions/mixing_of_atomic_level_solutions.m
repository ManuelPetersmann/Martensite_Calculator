function block_solutions = mixing_of_atomic_level_solutions(lath_solutions, block_solutions) % Solution_array_composite 
% call: mixing_of_atomic_level_solutions(atomic_solutions, tolerances)
% This function takes deformations F_i of atomic habit plane solutions and averages
% solutions fullfilling compatibility critera given by:

% tolerances... array of tolerance angles allowing for block mixing
% permitting the calculation of an average - following a linear rule of mixture: 
% prop_composite = x * prop_lath_sol1 + (1-x) prop_lath_sol2     e.g.
% F_composite = x*F_lath1 + (1-x)*F_lath2
% with x being the volume fraction of solution 1.
% Note, that for such an average the IPS condition is approximately
% maintained, which has been checked to be true for various combinations.
%
% The following cost-functions are minimized w.r.t the phase fraction x, 1-x 
% and are stored in an 'Composite_solution' object togehter with the invariant plane
% variables of the averaged block and the lath solutions used.
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


calculation_method = 'NEW Approach: Build blocks from lath-IPS-solutions, optimized phase fractions';
block_solutions.calculation_method = calculation_method;

I = eye(3); % = austenite
isol = 0;
neglected_due_to_eigenvalue_change_during_mixing = 0;
% loop over slip system combinations
for is1 = 1:(size(lath_solutions.array,2)-1)
    for is2 = (is1+1):size(lath_solutions.array,2)
        sol1 = lath_solutions.array(is1);
        sol2 = lath_solutions.array(is2);
        
        % do not mix variants not fullfilling predefined criteria
        % e.g. habit plane deviation from {111}_aust or something else
        % up to now there are two criteria ( both angle tolerances )
        % first angle between common line of invariant habit planes and
        % preferred invariant line (e.g. of set of cp-directions)
        % second angle between habit planes
        if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
            % considering longitudinal dimension of lath -> a
            vec_in_both_planes = cross( sol1.h , sol2.h );
            theta_intersec_cpdir = min_misorientation( lath_solutions.cryst_fams('KS'), vec_in_both_planes );
            if block_solutions.mixing_tolerances('theta_intersec_cpdir') < theta_intersec_cpdir
                continue
            end
        end
        %
        if isKey(block_solutions.mixing_tolerances,'theta_hps')
            theta_hps = get_angle( sol1.h , sol2.h ); 
            if block_solutions.mixing_tolerances('theta_hps') < theta_hps
            % considering width of laths -> b
            continue
            end
        end

        %if strcmp(minimize_prop,' shape strain vector')
        [x_eps, eps_block] = mixture_vecs_lin_least_squares_opt( cat(2,sol1.eps_ips*sol1.d , sol2.eps_ips*sol2.d) );
        F_comp_eps = reshape( cat(2,reshape(sol1.F1,9,1),reshape(sol1.F1,9,1)) * x_eps , 3,3) ;
        % it has been found that only for F_comp_eps the determinant is
        % invariant! 
        %
        [x_dis, d_dis] = mixture_matrix_lin_least_squares_opt(cat(3,sol1.ST,sol2.ST));
        % removed third output of above fucntion - ,  F_comp_dis_opt] since
        % determinant is not invariant - same for function call below
        [x_gl,  d_gl] = frob_min_green_lagrange_composite_block( sol1.ST,sol2.ST, x_eps(1) );
         
        % calculate composite IPS properties of F_comp_eps
        try
            [y1, y3, d1, d2, h1, h2, Q1, Q2] = rank_one(F_comp_eps, I, 1.e-3 ); 
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
        catch
            neglected_due_to_eigenvalue_change_during_mixing = neglected_due_to_eigenvalue_change_during_mixing + 2;
            continue % do not safe solution if
        end
        
        isol = isol + 2;
        % Create Slip_solution objects and append them to object array  - %
        % here i put LT = zeros(3) because it is not directly calculable!
        block_solutions.array( isol-1 ) =  Composite_solution(F_comp_eps, I, y1, y3, d1, h1, Q1, zeros(3),...
        eps_block, x_eps, F_comp_eps, x_dis, d_dis, x_gl, d_gl );
        block_solutions.array( isol )   =  Composite_solution(F_comp_eps, I, y1, y3, d2, h2, Q2, zeros(3),...
        eps_block, x_eps, F_comp_eps, x_dis, d_dis, x_gl, d_gl );
        %
        block_solutions.array( isol-1 ).lath_id_pair = [sol1.id, sol2.id];
        block_solutions.array( isol ).lath_id_pair =   [sol1.id, sol2.id];
        %
        if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
            block_solutions.array( isol-1 ).tolerances('theta_intersec_cpdir')  = theta_intersec_cpdir;
            block_solutions.array( isol ).tolerances('theta_intersec_cpdir')    = theta_intersec_cpdir;
        end
        if isKey(block_solutions.mixing_tolerances,'theta_hps')
            block_solutions.array( isol-1 ).tolerances('theta_hps') = theta_hps;
            block_solutions.array( isol ).tolerances('theta_hps')   = theta_hps;
        end
        %
    end % end of loop 1
end % end of loop 2

disp( ['number of potential solutions found: n_sol = ', num2str(isol) ] )
disp( [ num2str(neglected_due_to_eigenvalue_change_during_mixing), ' solutions neglected due to change of middle eigenvalue > 1.e-3'])
end