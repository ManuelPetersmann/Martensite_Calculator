function block_solutions = mixing_of_atomic_level_solutions(lath_solutions, U, cof_tol, det_tol) 

% call: mixing_of_atomic_level_solutions(lath_solutions, block_solutions, U, cof_tol, det_tol) 
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
    % minors tolerances
    cof_tol = 1.e-3
    det_tol = 1.e-3
end

block_solutions = Solution_array_composite();

I = eye(3); % = austenite
detU = det(U);
block_sols = 0;

rot_angle_block = 1.
%lambda2_tol_block_aust = 1.e-3 % doesnt matter if 0.001 or 0.0001 !!! important! some more solutions with 0.003
%block_hp_cp_aust_tol = 5.; % degree - even if i just set this only to 10 most solutions fall out
%lambda2_tol_laths = 1.e-4

neg_minors = 0;
neg_rot_angle = 0;
%neg_lamda2_block_aust = 0;
%neg_lamda2_laths = 0;
%neg_hp = 0;

% loop over slip system combinations
for is1 = 1: (size(lath_solutions.array,2)-1)
    for is2 = (is1+1): size(lath_solutions.array,2)
        sol1 = lath_solutions.array(is1);
        sol2 = lath_solutions.array(is2);
        F1 = sol1.ST;
        F2 = sol2.ST;
        
        x = 0.5;
        Fc = linmix2(x,F1,F2);
        
        %% Minors rules
        % third MINORS RULE
        det_Fc = det( Fc ); % plotting showed that if the determinant changes then the maximum deviation is at xi=0.5
        if abs(detU - det_Fc) > det_tol
            neg_minors = neg_minors + 1;
            continue
        end
        
        %% second MINORS RULE
        cofFc = cofactor( Fc );
        cof_F_sum = x * cofactor(F1)  +  (1.-x) * cofactor(F2);
        if sum(sum(abs(cofFc - cof_F_sum))) > cof_tol % alternatively frob distance
            neg_minors = neg_minors + 1;
            continue
        end
        
                 
%%       %% rotation of Block_inclusion
        [~,R] = polardecomposition( Fc );
        abs_angle = acosd( (trace(R)-1.) / 2.); 
        if abs_angle > rot_angle_block
            neg_rot_angle = neg_rot_angle +1;
            continue
        end

        %% ID pairs are all that is needed to form blocks and get all other information of them 
 
        block_sols = block_sols + 1;            
        block_solutions.array( block_sols ).lath_solution_pair = [sol1, sol2];  % U,tolerance]; %
        block_solutions.array( block_sols ).Fc05 = Fc;
        
        %% Currently no further a-priori mixing restrictions beside the minors conditions are implemented
%         if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
%         %    block_solutions.array( block_sols ).tolerances('theta_intersec_cpdir')    = theta_intersec_cpdir;
%         end
%         if isKey(block_solutions.mixing_tolerances,'theta_hps')
%         %    block_solutions.array( block_sols ).tolerances('theta_hps')   = theta_hps;
%         end
        
        
    end % end of loop 1
end % end of loop 2

combinations
disp( ['First  crit: ', num2str(neg_minors), ' neglected due to minors relations (cof_tol =',num2str(cof_tol), ', det_tol =',num2str(det_tol),')'] );
neg_rot_angle
%disp( ['Second crit: ', num2str(neg_lamda2_block_aust), ' neglected because lamda2 of block deviates more than ', num2str(lambda2_tol_block_aust), ' from 1'] );
%disp( ['Third crit: ', num2str(neg_lamda2_laths), ' neglected laths deformations of pairings are not rank one connected with tolerance ', num2str(lambda2_tol_laths) ] ); 
%disp( ['Fourth crit: ', num2str(neg_hp), ' neglected because ave HP deviates more than ', num2str(block_hp_cp_aust_tol),' from 111_aust']);
%disp( ['Fourth crit: ', num2str(neg_diff), ' neglected because Fs do not differ in the d1 norm more than ', num2str(delta_F_min) ] ); 

disp( ['number of potential solutions found: n_sol = ', num2str(block_sols) ] )

end