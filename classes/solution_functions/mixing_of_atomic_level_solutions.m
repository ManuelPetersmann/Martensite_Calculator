function block_solutions = mixing_of_atomic_level_solutions(lath_solutions, block_solutions,tol) % now i optimize for everything simultaneously, opt_func)  
% call: mixing_of_atomic_level_solutions(lath_solutions, block_solutions, opt_func)  
%
% lath_solutions ... array of lath solutions for building blocks
% block_solutions ... object of Solution_array_composite a.o. with property
% mixing_tolerances... dict/hashtable of tolerance angles allowing for block mixing
% permitting the calculation of an average - following a linear rule of mixture: 
% prop_composite = x * prop_lath_sol1 + (1-x) prop_lath_sol2     e.g.
% F_composite = x*F_lath1 + (1-x)*F_lath2
% with x being the volume fraction of solution 1.
% Note, that for such an average the IPS condition is approximately
% maintained, which has been checked to be true for various combinations.
%
% The following cost-functions 'opt_func' could be minimized w.r.t the phase fraction x, 1-x 
% and are stored in an 'Composite_solution' object togehter with the invariant plane
% variables of the averaged block and the lath solutions used.
%
% 1 ) 'eps': the magnitude of the linearly mixed shape strain vector as a measure of
% self accommodation.  % F_composite = I + eps * d_composite \otimes h_composite
% min_norm2( eps1*d1 + eps2*d2)--> linear constrained optimization
%
% 2 ) 'disg':  | F_composite -I |                   - displacement gradient
% 3 ) 'gl':    | F_composite^T F_composite - I |    - c.f. Green-Lagrange (no rotation)
%
% -)F_composite should be mostly volumetric, i.e.  min|F_composite - % F^spherical| --> nonlinear constrained optimization
% where the diagonal entries of F^spherical are determined from the volume 
% change (det(Bain)) - used e.g. in Qi2014 - though I think this is not reasonable


% Note: THE AVERAGE OF TWO IPS ALWAYS IS AN IPS (except for the small 
% error that is made because the determinant is not invariant to addition)
% In general the linear mixture rule is only valid in the case of
% n->infinity (minor relations) see Bhattacharya - Microstructures of
% martensties - p.131.

if nargin < 3
    tol = 1.e-3
end

calculation_method = 'NEW Approach: Build blocks from lath-IPS-solutions, optimized phase fractions';
block_solutions.calculation_method = calculation_method;

I = eye(3); % = austenite
isol = 0;
neglected_eigenvalue_change_during_mixing_opt_eps = 0;
neglected_eigenvalue_change_during_mixing_opt_dis = 0;
neglected_eigenvalue_change_during_mixing_opt_gl = 0;
neglected_mixing_restrictions = 0;
% loop over slip system combinations
for is1 = 1: 5 %(size(lath_solutions.array,2)-1)
    for is2 = (is1+1): 5 %size(lath_solutions.array,2)
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
                neglected_mixing_restrictions = neglected_mixing_restrictions +1;
                continue
            end
        end
        %
        if isKey(block_solutions.mixing_tolerances,'theta_hps')
            theta_hps = get_angle( sol1.h , sol2.h );
            if block_solutions.mixing_tolerances('theta_hps') < theta_hps
                % considering width of laths -> b
                neglected_mixing_restrictions = neglected_mixing_restrictions +1;
                continue
            end
        end
        
        % x = fmincon(@myfun,x0,A,b)
        
        %[ sol1.id , sol2.id ]
        
        [x_eps, eps_block] = mixture_vecs_lin_least_squares_opt( cat(2,sol1.eps_ips*sol1.d , sol2.eps_ips*sol2.d) );
        F_comp_eps = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_eps , 3,3) ;
        %x_eps
        %disp(['det(F_comp_eps) = ',num2str( det(F_comp_eps) ) ] );
        
        [x_dis, d_dis] = mixture_matrix_lin_least_squares_opt(cat(3,sol1.ST,sol2.ST));
        F_comp_dis = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_dis , 3,3) ;
        %x_dis
        %disp(['det(F_comp_dis) = ',num2str( det(F_comp_dis) ) ] );
        
        [x_gl,  d_gl] = frob_min_green_lagrange_composite_block( sol1.ST,sol2.ST, x_eps(1) );
        F_comp_gl = reshape( cat(2,reshape(sol1.ST,9,1),reshape(sol2.ST,9,1)) * x_gl , 3,3) ;
        %x_gl
        %disp(['det(F_comp_gl) = ',num2str( det(F_comp_dis) ) ] );
        
        %         switch opt_func
        %             case 'eps'
        %                 F_comp = F_comp_eps;
        %             case 'disg'
        %                 F_comp = F_comp_dis;
        %             case 'gl'
        %                 F_comp = F_comp_gl;
        %         end
        
        no_sol = {0, 0, zeros(1,3), zeros(1,3), zeros(1,3), zeros(1,3), zeros(3), zeros(3) };
        
        y1 = [ y1_eps;  y1_dis;  y1_gl ], ...
        y3 = [ y3_eps;  y3_dis;  y3_gl ], ...
        d = [ d1_eps;  d2_eps;  d1_dis;  d2_dis; d1_gl;  d2_gl ], ...
        h = [ h1_eps;  h2_eps;  h1_dis;  h2_dis; h1_gl;  h2_gl ], cat(3,Q1_eps, Q2_eps, Q1_dis, Q2_dis, Q1_gl, Q2_gl)
        
        % calculate composite IPS properties of all possibilites of optimized mixings
        try
            [y1_eps, y3_eps, d1_eps, d2_eps, h1_eps, h2_eps, Q1_eps, Q2_eps] = rank_one(F_comp_eps, I, tol); 
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
            
        catch %ME
            neglected_eigenvalue_change_during_mixing_opt_eps = neglected_eigenvalue_change_during_mixing_opt_eps + 2;
            %disp( ME );
            %[y1_eps, y3_eps, d1_eps, d2_eps, h1_eps, h2_eps, Q1_eps, Q2_eps] = deal(no_sol{:});
            %continue % do not safe solution if
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            [y1_dis, y3_dis, d1_dis, d2_dis, h1_dis, h2_dis, Q1_dis, Q2_dis] = rank_one(F_comp_dis, I, tol); 
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
        catch %ME
            neglected_eigenvalue_change_during_mixing_opt_dis = neglected_eigenvalue_change_during_mixing_opt_dis + 2;
            %disp( ME );
            %[y1_dis, y3_dis, d1_dis, d2_dis, h1_dis, h2_dis, Q1_dis, Q2_dis] = deal(no_sol{:});
            %continue % do not safe solution if
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            [y1_gl, y3_gl, d1_gl, d2_gl, h1_gl, h2_gl, Q1_gl, Q2_gl] = rank_one(F_comp_gl, I, tol); 
            % here a smaller tolerance value is used because the eigenvalue
            % changes slightly see introductory Note. If even this
            % tolerance is not reached the mixing is neglected
        catch %ME
            neglected_eigenvalue_change_during_mixing_opt_gl = neglected_eigenvalue_change_during_mixing_opt_gl + 2;
            %disp( ME );
            %[y1_gl, y3_gl, d1_gl, d2_gl, h1_gl, h2_gl, Q1_gl, Q2_gl] = deal(no_sol{:});
            %continue % do not safe solution if
        end
        
        isol = isol + 1;
        % Create Slip_solution objects and append them to object array
        % here i put LT = zeros(3) because it is not directly calculable!
        block_solutions.array( isol )   =  Composite_solution(F1, I, y1, y3, d, h, Q, ...
           zeros(3), eps_block, x_eps, F1, x_dis, d_dis, x_gl, d_gl );
            
%             [ y1_eps;  y1_dis;  y1_gl ], ...
%             [ y3_eps;  y3_dis;  y3_gl ],...
%             [ d1_eps;  d2_eps;  d1_dis;  d2_dis; d1_gl;  d2_gl ], ...
%             [ h1_eps;  h2_eps;  h1_dis;  h2_dis; h1_gl;  h2_gl ], cat(3,Q1_eps, Q2_eps, Q1_dis, Q2_dis, Q1_gl, Q2_gl), ...
            
        %
        block_solutions.array( isol ).lath_id_pair =   [sol1.id, sol2.id];
        %
        if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
            block_solutions.array( isol ).tolerances('theta_intersec_cpdir')    = theta_intersec_cpdir;
        end
        if isKey(block_solutions.mixing_tolerances,'theta_hps')
            block_solutions.array( isol ).tolerances('theta_hps')   = theta_hps;
        end
        %
    end % end of loop 1
end % end of loop 2

disp( [ num2str(neglected_mixing_restrictions), ' solutions selected out due restriction on mixing laths to blocks'])
disp( [ num2str(neglected_eigenvalue_change_during_mixing_opt_eps), ' solutions neglected due to change of middle eigenvalue > 1.e-3'] );
disp( [ num2str(neglected_eigenvalue_change_during_mixing_opt_dis),' solutions neglected due to change of middle eigenvalue > 1.e-3'] );
disp( [ num2str(neglected_eigenvalue_change_during_mixing_opt_gl) ,' solutions neglected due to change of middle eigenvalue > 1.e-3'] );
disp( ['number of potential solutions found: n_sol = ', num2str(isol) ] )


end