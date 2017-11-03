function block_solutions = mixing_of_atomic_level_solutions(lath_solutions, block_solutions, U, lambda2_tol, cof_tol, det_tol) % outarg - block_solutions tol) 
% call: mixing_of_atomic_level_solutions(lath_solutions, block_solutions, tol)  % opt_func  
%
% lath_solutions ... array of lath solutions for building blocks
% block_solutions ... object of Solution_array_composite a.o. with property
%
% block_solutions.mixing_tolerances... dict/hashtable of tolerance angles allowing for block mixing
%
% calculation of an average - following a linear rule of mixture: 
% prop_composite = x * prop_lath_sol1 + (1-x) prop_lath_sol2     e.g.
% F_composite = x*F_lath1 + (1-x)*F_lath2
% with x being the volume fraction of solution 1.
% Note, that for such an average the IPS condition is approximately
% maintained, which has been checked to be true for various combinations.
%
% Note: THE AVERAGE OF TWO IPS ALWAYS IS AN IPS (except for the small 
% error that is made because the determinant is not invariant to addition)
% In general the linear mixture rule is only valid in the case of
% n->infinity (minor relations) see Bhattacharya - Microstructures of
% martensties - p.131.

% now i optimize for everything simultaneously, opt_func)  

if nargin < 4
    lambda2_tol = 1.e-5;
    cof_tol = 1.e-5;
    det_tol = 1.e-5;
end

    function lambda2_mix = mix_y2( x, F1, F2)
        Fc = linmix2(x, F1, F2);
        [~,lambda2_mix] = sorted_eig_vals_and_vecs(Fc'*Fc);
    end

calculation_method = 'NEW Approach: Build blocks from lath-IPS-solutions, optimized phase fractions';
block_solutions.calculation_method = calculation_method;

I = eye(3); % = austenite
detU = det(U);
block_sols = 0;

% loop over slip system combinations
for is1 = 1: (size(lath_solutions.array,2)-1)
    for is2 = (is1+1): size(lath_solutions.array,2)
        sol1 = lath_solutions.array(is1);
        sol2 = lath_solutions.array(is2);
        F1 = sol1.ST;
        F2 = sol2.ST;
        
        x = 0.5;
        Fc = linmix2(x,F1,F2);
        
        %% third MINORS RULE
        det_Fc = det( Fc ); % plotting showed that if the determinant changes then the maximum deviation is at xi=0.5
        if abs(detU - det_Fc) > det_tol
            continue
        end
        
        %% second MINORS RULE
        cofFc = cofactor( Fc );
        cof_F_sum = x * cofactor(F1)  +  (1.-x) * cofactor(F2);
        if sum(sum(abs(cofFc - cof_F_sum))) > cof_tol % alternatively frob distance
            continue
        end
        
        %% check deviation of lambda2
        if (mix_y2(x,F1,F2) - 1)  >  lambda2_tol  
            continue
        end
        
        %% Considering that the two F could be equal since the slip
        % deformations are not linearly independent! c.f. non-uniqueness of
        % plastic slip
        if sum(sum(abs(F1 - F2 ))) < delta_F % alternatively frob distance
            continue
        end
        
        
        % do not mix variants not fullfilling predefined criteria
        % e.g. habit plane deviation from {111}_aust or something else
        % up to now there are two criteria ( both angle tolerances )
        % first angle between common line of invariant habit planes and
        % preferred invariant line (e.g. of set of cp-directions)
        % second angle between habit planes
        if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
            % considering longitudinal dimension of lath -> a
            vec_in_both_planes = cross( sol1.h , sol2.h );
            % check if habit planes are not parallel
            if abs(vec_in_both_planes) < 1.e-8 % entry wise for all entries
                theta_intersec_cpdir = misorientation_vector_and_plane( lath_solutions.cryst_fams('KS'), sol1.h );
            else
                theta_intersec_cpdir = min_misorientation( lath_solutions.cryst_fams('KS'), vec_in_both_planes );
            end
            %
            if theta_intersec_cpdir  >  block_solutions.mixing_tolerances('theta_intersec_cpdir')
                neglected_mixing_restrictions = neglected_mixing_restrictions +1;
                continue
            end
        end                 
        %
        if isKey(block_solutions.mixing_tolerances,'theta_hps')
            theta_hps = get_angle( sol1.h , sol2.h );
            if theta_hps  >  block_solutions.mixing_tolerances('theta_hps') 
                % considering width of laths -> b
                neglected_mixing_restrictions = neglected_mixing_restrictions +1;
                continue
            end
        end
        
        %% ID pairs are all that is needed to form blocks and get all other information of them 
 
        block_sols = block_sols + 1;            

        block_solutions.array( block_sols ).lath_id_pair =   [sol1.id, sol2.id];
        %
        if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
        %    block_solutions.array( block_sols ).tolerances('theta_intersec_cpdir')    = theta_intersec_cpdir;
        end
        if isKey(block_solutions.mixing_tolerances,'theta_hps')
        %    block_solutions.array( block_sols ).tolerances('theta_hps')   = theta_hps;
        end
        
        
    end % end of loop 1
end % end of loop 2


disp( ['number of potential solutions found: n_sol = ', num2str(block_sols) ] )


end