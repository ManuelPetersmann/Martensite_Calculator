function [min_sols, blocks] = minors_relation_and_IPS(lath_solutions, block_solutions, lambda2_tol, cof_tol, det_tol) % outarg - block_solutions
% call: minors_relation_and_IPS(lath_solutions, block_solutions,tol)
%
% lath_solutions ... array of lath solutions for building blocks
% block_solutions ... object of Solution_array_composite a.o. with property

% F_composite = x*F_lath1 + (1-x)*F_lath2
% with x being the volume fraction of solution 1.
% Note, that for such an average the IPS condition is approximately
% maintained, which has been checked to be true for various combinations.

% Note: THE AVERAGE OF TWO IPS ALWAYS IS AN IPS (except for the small 
% error that is made because the determinant is not invariant to addition)
% In general the linear mixture rule is only valid in the case of
% n->infinity (minor relations) see Bhattacharya - Microstructures of
% martensties - p.131.

if nargin < 3
    lambda2_tol = 1.e-8;
    cof_tol = 1.e-5;
    det_tol = 1.e-5;
end

    function y2 = mix_y2(x,F1,F2)
        Fc = x * F1  +  (1.-x) * F2;
        [ ~, y2 ] = sorted_eig_vals_and_vecs( Fc'*Fc );
    end

%calculation_method = 'Minors relations two matrices';
%block_solutions.calculation_method = calculation_method;

I = eye(3); % = austenite
isol = 1;
min_sols = 0;
blocks = 0;

% neglected_eigenvalue_change_during_mixing_opt_eps = 0;
% neglected_eigenvalue_change_during_mixing_opt_dis = 0;
% neglected_eigenvalue_change_during_mixing_opt_gl = 0;
% neglected_mixing_restrictions = 0;

%xi = linspace(0,1,100);
xi = linspace(0.1,0.85,6);
%xi = linspace(0,1,2);
y2 = zeros(1,20);

% loop over slip system combinations
for is1 = 10: 20 %(size(lath_solutions.array,2)-1)
    for is2 = 20 %(is1+1):size(lath_solutions.array,2)
        sol1 = lath_solutions.array(is1);
        sol2 = lath_solutions.array(is2);
        
        % do not mix variants not fullfilling predefined criteria
        % e.g. habit plane deviation from {111}_aust or something else
        % up to now there are two criteria ( both angle tolerances )
        % first angle between common line of invariant habit planes and
        % preferred invariant line (e.g. of set of cp-directions)
        % second angle between habit planes
        %         if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
        %             % considering longitudinal dimension of lath -> a
        %             vec_in_both_planes = cross( sol1.h , sol2.h );
        %             theta_intersec_cpdir = min_misorientation( lath_solutions.cryst_fams('KS'), vec_in_both_planes );
        %             if block_solutions.mixing_tolerances('theta_intersec_cpdir') < theta_intersec_cpdir
        %                 neglected_mixing_restrictions = neglected_mixing_restrictions +1;
        %                 continue
        %             end
        %         end
        %         %
        %         if isKey(block_solutions.mixing_tolerances,'theta_hps')
        %             theta_hps = get_angle( sol1.h , sol2.h );
        %             if block_solutions.mixing_tolerances('theta_hps') < theta_hps
        %                 % considering width of laths -> b
        %                 neglected_mixing_restrictions = neglected_mixing_restrictions +1;
        %                 continue
        %             end
        %         end

        F1 = sol1.ST;
        F2 = sol2.ST; 
        
        % reformulation of FIRST MINORS RULE - IPS constraint on the boundary:
        % Bisection to find solution for xi between 0 and 1.
        x_left  = 0.02; % linke intervallschranke -> 1-x rechte intervallschranke
        x_right = 1.- x_left; 
        y2_left  = mix_y2( x_left, F1, F2);
        y2_right = mix_y2( x_right, F1, F2);
        % yy = [y2_1, y2_2];   
        %if ( any( yy > 1)  &&  any( yy < 1) ) 
        if (y2_right -1. ) * (y2_left -1. ) < 0
            while true
                x = (x_right + x_left)/2.;
                mid = mix_y2( x, F1, F2) ;
                if abs(mid -1.) < lambda2_tol
                %    fracs(isol) = x;
                    blocks = blocks +1;
                    break
                end
                %
                if (y2_right -1. ) * (mid -1. ) < 0
                    x_left = x;
                else
                    x_right = x;
                    y2_right = mix_y2( x_right, F1, F2);
                end
            end
        else % no local minimum for block ! Do not mix laths to blocks
            continue
        end
        
        if mod(blocks,1000)==0
            blocks
        end
        
        xi(1) = x;
        for i=1:length(xi) 
            y2(i) = mix_y2( xi(i), F1, F2);
            figure;
            plot(xi,y2);
            x = xi(i);
            
%             % second MINORS RULE
%             cof_F_ips = cofactor( x * F1  +  (1.-x) * F2 );
%             cof_F_sum = x * cofactor(F1)  +  (1.-x) * cofactor(F2);
%             if frob_distance(cof_F_ips , cof_F_sum) > cof_tol
%                 continue
%             end
%             
%             % third MINORS RULE
%             det_F_ips = det( x * F1  +  (1.-x) * F2 );
%             det_F_sum = x * det(F1)  +  (1.-x) * det(F2);
%             if abs(det_F_ips - det_F_sum) > det_tol
%                 continue
%             end
%             
%             if i == 1
%                 disp('composite frac check')
%                 xi(1)
%             end
%             x
%             min_sols = min_sols + 1;
            
        end  % end for
        

        
%         while true
%             % calculate composite IPS properties of all possibilites of optimized mixings
%             try
%                 [y1_eps, y3_eps, d1_eps, d2_eps, h1_eps, h2_eps, Q1_eps, Q2_eps] = rank_one(F_comp_eps, I, tol);
%                 % here a smaller tolerance value is used because the eigenvalue
%                 % changes slightly see introductory Note. If even this
%                 % tolerance is not reached the mixing is neglected
%             catch %ME
%                 neglected_eigenvalue_change_during_mixing_opt_eps = neglected_eigenvalue_change_during_mixing_opt_eps + 2;
%                 %disp( ME );
%                 %[y1_eps, y3_eps, d1_eps, d2_eps, h1_eps, h2_eps, Q1_eps, Q2_eps] = deal(no_sol{:});
%                 %continue % do not safe solution if
%             end
%             
%         end
        
%         isol = isol + 1;
%         % Create Slip_solution objects and append them to object array
%         % here i put LT = zeros(3) because it is not directly calculable!
%         block_solutions.array( isol )   =  Composite_solution(F1, I, y1, y3, d, h, Q, ...
%             zeros(3), eps_block, x_eps, F1, x_dis, d_dis, x_gl, d_gl );
%         %
%         block_solutions.array( isol ).lath_id_pair =   [sol1.id, sol2.id];
        %
        %         if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
        %             block_solutions.array( isol ).tolerances('theta_intersec_cpdir')    = theta_intersec_cpdir;
        %         end
        %         if isKey(block_solutions.mixing_tolerances,'theta_hps')
        %             block_solutions.array( isol ).tolerances('theta_hps')   = theta_hps;
        %         end
        %
    end % end of loop 1
end % end of loop 2

% disp( [ num2str(neglected_mixing_restrictions), ' solutions selected out due restriction on mixing laths to blocks'])
% disp( [ num2str(neglected_eigenvalue_change_during_mixing_opt_eps), ' solutions neglected due to change of middle eigenvalue > 1.e-3'] );
% disp( [ num2str(neglected_eigenvalue_change_during_mixing_opt_dis),' solutions neglected due to change of middle eigenvalue > 1.e-3'] );
% disp( [ num2str(neglected_eigenvalue_change_during_mixing_opt_gl) ,' solutions neglected due to change of middle eigenvalue > 1.e-3'] );

disp( ['number of potential solutions found: n_sol = ', num2str(isol) ] )


end