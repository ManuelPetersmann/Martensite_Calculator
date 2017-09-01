function [solutions] = doubleshear_variable_shear_mags(B, ns_product, ds_product, cp, ns_parent, ds_parent)
% possible calls: multiple_shears_incremental_optimization(B, ns_parent, ds_parent)
%                                                         (B,ns_product, ds_product, cp)
%                                                         (B, ns_product, ds_product, cp, ns_parent, ds_parent) 
% Function can be called with 3 (only parent slip systems) 4 (only product slip systems), 6 (slip systems of both phases)
% All calulations are carried out in the coordinate system of the parent phase
% B... Bain strain, 
% cp... B*Correspondance matrix --- mapping parent phase vectors to product phase vectors
% ns...slip system normals, ds... slip directios
% returns object array of solutions for IPSs.

solutions = Solution_array( Slip_solution_doubleshear() ); 
% Construct array with type of solution -> After this line, Solution_array.array is no longer a double 

%% set numerical parameters (see file numerical_parameters.m)
numerical_parameters;

%% transform product phase slip systems to parent phase and combine all in one array
%ds = zeros(1,3);
%ns = zeros(1,3);
if nargin == 3 % only parent phase slip systems
    ds = ds_product;
    ns = ns_product;
end
if nargin > 3
    for is = 1:size(ds_product,1)
        % transform product phase slip systems to parent phase ones
        ds(is,:) = cp * ds_product(is,:)';
        ns(is,:) = inverse(cp)' * ns_product(is,:)';
    end
end
if nargin == 6  % if both parent and product phase systems are given
    ds = cat(1,ds,ds_parent);
    ns = cat(1,ns,ns_parent);
end

% no normed vectors in Khachaturyans approach!!
for jj = 1:size(ds,1)
%    ds(jj,:) = ds(jj,:) / norm(ds(jj,:));
%    ns(jj,:) = ns(jj,:) / norm(ns(jj,:));
    S(:,:,jj)  = (ds(jj,:)' * ns(jj,:) ); 
end

%% calculate only initial eigenvalues without shear modification to determine
% the direction from which side lambda2 = 1 is approached
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( B'*B );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, epsilon);
old_min_delta_lambda2_to_1 = abs(1. - lambda_2);

% loop over slip system combinations
for is1 = 1:(size(ds,1)-1) % loop for first slip system
    for is2 = (is1+1):size(ds,1) % loop for second one              
        % modify shear value in Blocks until lambda2 = 1
        %ds(is2,:)
        %ns(is2,:)
        delta_g = delta_g_initial;
        g1 = g_initial;
        g2 = g_initial;
        % g = g_initial / 1./(norm(d11)*norm(n11));
        is_possible_solution = false;
        %lambda2_smaller1 = lambda2_smaller1_initial;
        while ( ~is_possible_solution && ((g1 > g_min) || g2 > g_min) )  % if the solution for g is very high or low respectively, do not consider it
            
            if ((g1 > g_initial) || (g2 > g_initial))
                error('this should not happen - fix code...')
            end
            
            S1 =  I + (1./(g1 - delta_g)) * S(:,:,is1);  % shear amplitude gets bigger with smaller g!
            F = B * S1;
            [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
            [ is_possible_solution , lambda2_smaller1_shear1 ] = check_IPS_solution(lambda_1, lambda_2, lambda_3, epsilon);
            if is_possible_solution
                break
            end
            delta_lambda2_S1 = abs(1. - lambda_2);
            
            S2 =  I + (1./(g2 - delta_g)) * S(:,:,is2);  % shear amplitude gets bigger with smaller g!
            F = B * S2; % here it has been tested that the order of multiplication does not matter
                        % since the Bain is pure stretch and SS is a small strain
            [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
            [ is_possible_solution , lambda2_smaller1_shear2] = check_IPS_solution(lambda_1, lambda_2, lambda_3, epsilon);
            if is_possible_solution
                break
            end
            delta_lambda2_S2 = abs(1. - lambda_2);
            
            
            % choose the shear that approached a solution quicker with the
            % SAME shear magnitude g
            if(delta_lambda2_S1 < delta_lambda2_S2) % here generally the minimum should be taken if the solution has not been passed already
                % from one search direction, otherwise the shear amplitude is adopted first and only then it is
                % checked again wheter a new minimum does not overshoot the solution from one direction
                if ( ( lambda2_smaller1_shear1 ~= lambda2_smaller1_initial ) || ( old_min_delta_lambda2_to_1 < delta_lambda2_S1 ) )
                    % i.e. the solution is passed or delta_lambda2_to_1 does not decrease, which should not happen - cut back
                    delta_g = stepwidth_change * delta_g; 
                else
                    old_min_delta_lambda2_to_1 = delta_lambda2_S1
                    g1 = g1 - delta_g; % shear amplitude gets bigger with smaller g!
                end
            else % (delta_lambda2_S1 > delta_lambda2_S2)
                if ( ( lambda2_smaller1_shear2 ~= lambda2_smaller1_initial ) || ( old_min_delta_lambda2_to_1 < delta_lambda2_S1 ) )
                    delta_g = stepwidth_change * delta_g;
                else
                    old_min_delta_lambda2_to_1 = delta_lambda2_S2
                    g2 = g2 - delta_g;
                end
            end
            
            if delta_g < 1.e-20
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
            [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F, I );
            % Note habit plane solutions come in pairs!
            
            isol = isol + 2 % increase counter for number of solutions found
            
            % Create Slip_solution objects and append them to object array
            solutions.array( isol-1 ) =  Slip_solution_doubleshear(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, g1, ds(is1,:), ns(is1,:), g2, ds(is2,:), ns(is2,:) );
            solutions.array( isol )   =  Slip_solution_doubleshear(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, g1, ds(is1,:), ns(is1,:), g2, ds(is2,:), ns(is2,:) );
        end
        
    end % end of loop for second slip system
end % end of loop for first slip system


fprintf('number of potential solutions found: n_sol = %i :\n', isol)

end


