function [solutions] = doubleshear(B, ns, ds, cp)
% call: oubleshear(B, cp, ms, ns, ds)
% B... Bain strain, cp - B*Correspondance matrix alpha',
% ns...slip system normals in alpha, ds... slip directios in alpha
% returns object array of solutions for IPSs.

if nargin < 4
    cp = eye(3); % if the slip systems are given 
end

%% initalize some other vars
solutions = Solution_array( Slip_solution() ); % Construct array with type of solution -> After this line, Solution_array.array is no longer a double 
epsilon = 1.e-9; % accuracy for middle valued eigenvalue
g_min = 4.;
g_initial = 100.0; 
delta_g_initial = 5.; % must be set reasonably so that a solution is found i.e. g_min not to high, delta_g not to high either
% here it reaches as lowest g=5 with delta_g = 5. everything lower will be disregarded
I = eye(3);
isol = 0; % counter for number of solutions

%% calculate only initial eigenvalues without shear modification to determine
% the direction from which side lambda2 = 1 is approached
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( B'*B );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, epsilon);

% loop over slip system combinations
for is1 = 1:(size(ds,1)-1) % loop for first slip system
    for is2 = (is1+1):size(ds,1) % loop for second one
        
        % ############  transform bcc systems to austenite ############
        % first system
        d1 = cp * ds(is1,:)'; % martensite.vec_from_coords( ds(is1,:) )'; - not necessary for cubic lattice
        n1 = inverse(cp)' * ns(is1,:)';

        % second system
        d2 = cp * ds(is2,:)';
        n2 = inverse(cp)' * ns(is2,:)';

        S1  = (d1  * n1') ;
        S2  = (d2  * n2') ;

        
        %% modify shear value in Blocks until lambda2 = 1
        delta_g = delta_g_initial;
        g = g_initial;
        % g = g_initial / 1./(norm(d11)*norm(n11));
        is_possible_solution = false;
        lambda2_smaller1 = lambda2_smaller1_initial;
        while ( ~is_possible_solution && (g > g_min) )  % if the solution for g is very high or low respectively, do not consider it
            
            if g > g_initial
                error('this should not happen - fix code...')
            end
            
            S =  I + (1./g)* (S1 + S2);

            % Construct the net deformation gradient from the two
            % sheared sides
            F = S * B; % composite block deformations, matrix multiplication is distributive
            % get new results
            [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
            
            %% check if solution has been found or how it changed if its not sufficient
            [ is_possible_solution , lambda2_smaller1_new] = check_IPS_solution(lambda_1, lambda_2, lambda_3, epsilon);
            
            % change the search direction and reduce step intervall if lambda2
            % passes one but is not in the required precision range.
            if lambda2_smaller1 ~= lambda2_smaller1_new
                delta_g = (-1./2.) * delta_g;
                %error('passed 1...')
            end
            lambda2_smaller1 = lambda2_smaller1_new;
            
            % change g value.
            g = g - delta_g;
            % find g (shear magnitude - m in Paper Qi, Khachaturyan 2014)
            % within specified limits (g_min, g_max)
            % e.g. g = 10.0 : 0.1 : 50.0  % maximum shear all 10 layers a step minmum all 50
            % negative g values are not necessary since mirror symmetry
            % is assumed and the solutions are already entailed
            
        end % end while
        
        
        if is_possible_solution
            %% calculate solution
            % calculate invariant plane vector n_i etc.
            [eps_0, a1, a2, n1, n2, Q1, Q2] = rank_one(F, I );
            % Note habit plane solutions come in pairs!
            
            isol = isol + 2; % increase counter for number of solutions found
            
            % Create Slip_solution objects and append them to object array
            solutions.array( isol-1 ) =  Slip_solution(F, I, isol-1, eps_0, a1, n1, Q1, Q1*B, g, ds(is1,:), ns(is1,:), ds(is2,:), ns(is2,:) );
            solutions.array( isol )   =  Slip_solution(F, I, isol,   eps_0, a2, n2, Q2, Q2*B, g, ds(is1,:), ns(is1,:), ds(is2,:), ns(is2,:) );
        end
        
    end % end of loop for second slip system
end % end of loop for first slip system



fprintf('number of potential solutions found: n_sol = %i :\n', isol)

end


