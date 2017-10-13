function doubleshear_variable_shear_mags(martensite, austenite) % [solutions] = 
% incremental optimization of "distance to middle eigenvalue = 1" approach on the lath level of highly dislocated lath martesite
% after Petersmann et al 2017 -Blocky lath martensite - theory,experiments and modeling - to be submitted
% All calulations are carried out in the coordinate system of the parent phase
% returns object array of solutions for IPSs.

% specify type of solution array
martensite.IPS_solutions.array = Slip_solution();
% set calcuation method property in solution_array object
calculation_method = 'variable doubleshear incremental optimization lath level';
martensite.IPS_solutions.calculation_method = calculation_method;  
% create shorthand notation
solutions = martensite.IPS_solutions;

%% set numerical parameters (see file numerical_parameters.m)
numerical_parameters;

%% transform product phase slip systems to parent phase and combine all in one array
% assemble all shear dyads in austenite, array of directions, planes and in
% respective phase (miller indizes)
[ds, ns, S, slip_combinations] = shear_dyads(martensite, austenite, false); % assemble normed- shear_dyads
solutions.slip_combinations = slip_combinations; % nr of possibilites nchoosek (k=2)

%% calculate only initial eigenvalues without shear modification to determine
% the direction from which side lambda2 = 1 is approached
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( martensite.U' * martensite.U );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, tolerance);
delta_lambda2_to_1_initial = abs(1. - lambda_2);

% loop over slip system combinations
for is1 = 1:(size(ds,1)-1) % loop for first slip system
    for is2 = (is1+1):size(ds,1) % loop for second one              
        % modify shear value in Blocks until lambda2 = 1
        %ds(is2,:)
        %ns(is2,:)
        old_min_delta_lambda2_to_1 = delta_lambda2_to_1_initial;
        delta_eps = delta_eps_initial;
        eps1 = eps_initial;
        eps2 = eps_initial;
        % g = g_initial / 1./(norm(d11)*norm(n11));
        is_possible_solution = false;
        S_accummulated = I;
        %lambda2_smaller1 = lambda2_smaller1_initial;
        while ( ~is_possible_solution && (eps1 < eps_max) && (eps2 < eps_max) )  
            % if the solution for g is very high or low respectively, do not consider it
            if ((eps1 < eps_initial) || (eps2 < eps_initial))
                error('this should not happen - fix code...')
            end
            
            S1 =  S_accummulated * (I + delta_eps* S(:,:,is1) );  
            F = martensite.U * S1;
            [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
            [ is_possible_solution , lambda2_smaller1_shear1 ] = check_IPS_solution(lambda_1, lambda_2, lambda_3, tolerance);
            if is_possible_solution
                break
            end
            new_delta_lambda2_S1 = abs(1. - lambda_2);
            
            S2 =  S_accummulated * (I + delta_eps* S(:,:,is2) );  
            F = martensite.U * S2; % here it has been tested that the order of multiplication does not matter
                        % since the Bain is pure stretch and SS is a small strain
            [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
            [ is_possible_solution , lambda2_smaller1_shear2] = check_IPS_solution(lambda_1, lambda_2, lambda_3, tolerance);
            if is_possible_solution
                break
            end
            new_delta_lambda2_S2 = abs(1. - lambda_2);           
            
            % choose the shear that approached a solution quicker with the
            % SAME shear magnitude g
            if(new_delta_lambda2_S1 < new_delta_lambda2_S2) 
                % here generally the minimum should be taken if the solution has not been passed already
                % from one search direction, otherwise the shear amplitude is adopted first and only then it is
                % checked again wheter a new minimum does not overshoot the solution from one direction
                if ( ( lambda2_smaller1_shear1 ~= lambda2_smaller1_initial ) || ( old_min_delta_lambda2_to_1 < new_delta_lambda2_S1 ) )
                    % i.e. the solution is passed or delta_lambda2_to_1 does not decrease, which should not happen - cut back
                    delta_eps = stepwidth_change * delta_eps;
                else
                    old_min_delta_lambda2_to_1 = new_delta_lambda2_S1;
                    eps1 = eps1 + delta_eps;
                    S_accummulated = S_accummulated * ( I + delta_eps * S(:,:,is1) );
                end
            else % (delta_lambda2_S1 > delta_lambda2_S2)
                if ( ( lambda2_smaller1_shear2 ~= lambda2_smaller1_initial ) || ( old_min_delta_lambda2_to_1 < new_delta_lambda2_S2 ) )
                    delta_eps = stepwidth_change * delta_eps;
                else
                    old_min_delta_lambda2_to_1 = new_delta_lambda2_S2;
                    eps2 = eps2 + delta_eps;
                    S_accummulated = S_accummulated * ( I + delta_eps * S(:,:,is2) );
                end
            end
            
            if delta_eps < tolerance 
                break
            end          
            % find g (shear magnitude - m in Paper Qi, Khachaturyan 2014)
            % within specified limits (g_min, g_max)
            % e.g. g = 10.0 : 0.1 : 50.0  % maximum shear all 10 layers a step minmum all 50
            % negative g values are not necessary since mirror symmetry
            % is assumed and the solutions are already entailed
            
        end % end while
              
        if is_possible_solution
            %% calculate solution
            % calculate invariant plane vector n_i etc.
            [y1,y3, d1, d2, h1, h2, Q1, Q2] = rank_one(F, I, tolerance );
            % Note habit plane solutions come in pairs!
            
            isol = isol + 2; % increase counter for number of solutions found
            if mod(isol,500)==0
                isol
                %pause(1);
            end
            eps_s = [eps1; eps2];
            d = [ds(is1,:); ds(is2,:)];
            n = [ns(is1,:); ns(is2,:)];
            
            % Create Slip_solution objects and append them to object array;
            % PET 10.10.17: replaced 'isol' and 'eps' wit y1 and y2
            solutions.array( isol-1 ) =  Slip_solution(F, I, y1, y3, d1, h1, Q1, Q1*martensite.U, eps_s, d, n );
            solutions.array( isol )   =  Slip_solution(F, I, y1, y3, d2, h2, Q2, Q2*martensite.U, eps_s, d ,n );
        end
        
    end % end of loop for second slip system
end % end of loop for first slip system


disp(['number of potential solutions found = ', num2str(isol)])

end


