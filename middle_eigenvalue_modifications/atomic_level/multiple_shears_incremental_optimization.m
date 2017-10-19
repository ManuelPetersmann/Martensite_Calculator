function [solutions] = multiple_shears_incremental_optimization(martensite, austenite) 
% call: multiple_shears_incremental_optimization(martensite, austenite)
% incrementally minimizes middle-valued eigenvalue over the whole set of
% selected slip systems
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
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( martensite.U'*martensite.U );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, tolerance);
old_min_delta_lambda2_to_1 = abs(1. - lambda_2);

%% modify shear value in Blocks incrementally until lambda2 = 1
shear_mag = zeros(1,size(ds,1)); % accummulated shears = 1/g
is_possible_solution = false;
S_accummulated = I;

% for printing evolution of certain shear magnitudes vs delta_lambda2-1
%delta_lambda2_to_1 = [];
%smag1 = [];
%smag2 = [];

lambda2_smaller1(1,size(ds,1)) = lambda2_smaller1_initial;
while ( ~is_possible_solution && ( all(shear_mag) < eps_max ) )
    % if the optimal solution yields a very low g value for
    % a specific system, do not consider it
    if ( any(shear_mag(:) > g_initial) )
        error('this should not happen - fix code...')
    end
    deltas_to_lambda2 = ones(1,size(ds,1))*5; % should not be larger than one but safety first
       
    for is = 1:size(ds,1)
        % add/test incremental/virtual shear and select the one that
        % converges the fastest against lambda_2 = 1
        SS =  S_accummulated * ( I + delta_eps * S(:,:,is) );  
        % here + or - plays no role because both shear directions are considered
        F = martensite.U*SS; % here it has been tested that the order of multiplication does not matter
                  % since the Bain is pure stretch and SS is a small strain
        [ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( F'*F );
        [ is_possible_solution , lambda2_smaller1(is) ] = check_IPS_solution(lambda_1, lambda_2, lambda_3, tolerance);
        if is_possible_solution
            break
        end
        deltas_to_lambda2(is) = abs(1. - lambda_2);
    end % end for
    if is_possible_solution
        break
    end
    
    % find minimum value
    [new_min_delta_lambda2_to_1, indexof_mind2l2] = min( deltas_to_lambda2 );
    new_min_delta_lambda2_to_1 %= min( deltas_to_lambda2 )
    
    % find all positons if there are more equal minima
    %new_min_delta_lambda2_to_1 =
    nr_of_minia = size( deltas_to_lambda2( abs(deltas_to_lambda2 - new_min_delta_lambda2_to_1) < 1.e-18 ) )
    % TODO: add as second incremental criterion the shape strain
    % lambda_1-lambda_3 to chose best if there are >1 nr_of_minima
    
    % choose the shear that approached a solution fastest in an increment (all with the SAME shear magnitude delta_eps)
    % here generally the minimum should be taken. If the solution has been passed already
    % from one search direction, the shear amplitude is adopted then it is
    % checked again wheter a new minimum does not overshoot the solution from one direction
    if ( any(lambda2_smaller1(:) ~= lambda2_smaller1_initial)  || ( old_min_delta_lambda2_to_1 < new_min_delta_lambda2_to_1 ) ) 
        % i.e. the solution is passed, which should not happen - cut back
        % or the distance to lambda_2=1 increased during increment
            delta_eps = stepwidth_change * delta_eps % 0.5 : make new incremental shear smaller!
    else
        shear_mag(indexof_mind2l2) =  shear_mag(indexof_mind2l2) + abs(delta_eps); % shear amplitude gets bigger with smaller g!
        % Udate S_accummulated
        S_accummulated = S_accummulated * ( I + delta_eps * S(:,:,indexof_mind2l2) );
        % throw error if local divergence occurs
        if abs(old_min_delta_lambda2_to_1 - new_min_delta_lambda2_to_1) < 1.e-18
            shear_mag % print out shear magnitudes before error is thrown
            error(['local divergence of middle_valued_eigenvalue at ' num2str(new_min_delta_lambda2_to_1), ' try more slipsystems or other numerical parameters'])
        else
            old_min_delta_lambda2_to_1 = new_min_delta_lambda2_to_1
            %
            % to print evolution of certain shear_magnitudes vs delta_lambda2_to_1
            %delta_lambda2_to_1( size(delta_lambda2_to_1,2) + 1) = old_min_delta_lambda2_to_1;
            %smag1 = []
            %smag2 = [];
        end
    end
    % find g (shear magnitude - m in Paper Qi, Khachaturyan 2014)
    % within specified limits (g_min, g_max)
    % e.g. g = 10.0 : 0.1 : 50.0  % maximum shear all 10 layers a step minmum all 50
    % negative g values are not necessary since each slip system is counted twice
     
end % end while


if is_possible_solution
    %% calculate solution
    % calculate invariant plane vector n_i etc.
    [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F, I, tolerance );
    % Note habit plane solutions come in pairs!
    
    isol = isol + 2; % increase counter for number of solutions found
    
    found = (shear_mag > 1.e-15);
    %shear_mag
    shear_magsfound = shear_mag(found);
    %
    dsfound = ds(found,:);
    nsfound = ns(found,:);
    
    % Create Slip_solution objects and append them to object array
    solutions.array( isol-1 ) =  Slip_solution_multishear(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, shear_magsfound, dsfound, nsfound );
    solutions.array( isol )   =  Slip_solution_multishear(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, shear_magsfound, dsfound, nsfound );
end

fprintf('number of potential solutions found: n_sol = %i :\n', isol)

end


