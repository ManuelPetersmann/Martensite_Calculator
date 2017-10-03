function [solutions] = block_symmetric_doubleshear(B, ms, ns_product, ds_product, cp, ns_parent, ds_parent)
% possible calls: block_symmetric_doubleshear(B, ms, ns_parent, ds_parent)
%                                                         (B, ms, ns_product, ds_product, cp)
%                                                         (B, ms, ns_product, ds_product, cp, ns_parent, ds_parent) 
% Function can be called with 3 (only parent slip systems) 4 (only product
% slip systems), 6 (slip systems of both phases) arguments
% All calulations are carried out in the coordinate system of the parent phase
% B... Bain strain, 
% cp... B*Correspondance matrix --- mapping parent phase vectors to product phase vectors
% ns...slip system normals, ds... slip directios
% ms ... mirror planes for construction of blocks from one Bain
% returns object array of solutions for IPSs.

numerical_parameters;
solutions = Solution_array( Slip_solution() );

%% transform product phase slip systems to parent phase and combine all in one array

% assemble all shear directions, planes and dyads
%[ds, ns, S] = shear_dyads(martensite.considered_plasticity, );

if nargin == 4 % only parent phase slip systems
    ds = ds_product;
    ns = ns_product;
end
if nargin > 4
    for is = 1:size(ds_product,1)
        % transform product phase slip systems to parent phase ones
        ds(is,:) = cp * ds_product(is,:)';
        ns(is,:) = inverse(cp)' * ns_product(is,:)';
    end
end
if nargin == 7  % if both parent and product phase systems are given
    ds = cat(1,ds,ds_parent);
    ns = cat(1,ns,ns_parent);
    % for outputting found slip systems in miller indizes
    ds_product = cat(1,ds_product,ds_parent);
    ns_product = cat(1,ns_product,ds_parent);
end

% to write integer values into solutions
for jj = 1:size(ds,1)
    % ds(jj,:) = ds(jj,:) / norm(ds(jj,:));
    % ns(jj,:) = ns(jj,:) / norm(ns(jj,:));
    % NO NORMATION!!! OTHER PREFACTOR !!!
    S(:,:,jj)  = ds(jj,:)' * ns(jj,:);
end

disp( ['Number of possible pairings is = ', num2str( nchoosek(size(ds,1),2) )])
disp('nr of solutions cannot be greater than 2-times this value.')


%% calculate only initial eigenvalues without shear modification to determine
% the direction from which side lambda2 = 1 is approached
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( B'*B );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, tolerance);

% lambda2_old = lambda_2;  % not used since it is not optimized
% incrementally...

%% loop over mirror planes and slip systems
for im = 1:size(ms,1) % number of considered mirror planes in martensite
    
    m_mart = ms(im,:); % mirror plane in martensite
    m_aust = inverse(cp)' * m_mart'; % transformed plane in austenite
    
    % loop over slip system combinations
    for is1 = 1:(size(ds,1)-1) % loop for first slip system
        for is2 = (is1+1):size(ds,1) % loop for second one
            
            d1_mirr = mirror_by_plane(m_aust, ds(is1,:), I);
            n1_mirr = mirror_by_plane(m_aust, ns(is1,:), I);
            d2_mirr = mirror_by_plane(m_aust, d2, I);
            n2_mirr = mirror_by_plane(m_aust, n2, I);
            S11 = (d1_mirr * n1_mirr') ;
            S22 = (d2_mirr * n2_mirr') ;
            
            %% modify shear value in Blocks until lambda2 = 1
            delta_g = delta_g_initial;
            g = g_initial;
            % g = g_initial / 1./(norm(d11)*norm(n11));
            is_possible_solution = false;
            lambda2_smaller1 = lambda2_smaller1_initial;
            while ( ~ is_possible_solution && (g > g_min) )  
                % if the solution for g is very high or low respectively,
                % do not consider it                
                if g > g_initial
                    error('this should not happen - fix code...')
                end
                
                % given the slip direction d and the slip plane normal n and the shear magnitude
                % 1/g (i.e. the lattice has a step after each g planes after the shear)
                % calculate double shear matrizes
                
                % shears are not commutative in large/finite strain: Sx = (I + S1)*(I + S2) \uneq Sy = (I + S2)*(I + S1)
                % Using a small/infinite strain assumption (reasonable
                % since slip should be small) the order does not matter.
                % S = (I + (1./g)*S1)*(I + (1./g)*S2)
                % S_mirror = (I + (1./g)*S11)*(I + (1./g)*S22)
                % verified that like Khachaturyan writes it, the order does not matter
                % i.e. the small strain assumption is justified!  
                
                S =  I + (1./g)* (S(:,:,is1) + S(:,:,is2));
                S_mirror = I + (1./g)* (S11 + S22);
                
                % calculate the rotation of the mirror plane vector due to the shear
                R = max_shear_rotation( m_aust, S);
                
                % Construct the net deformation gradient from the two
                % sheared sides
                F = 0.5*( R*S + inverse(R)*S_mirror ) * B; % composite block deformations, matrix multiplication is distributive
                
                % get new results
                [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
                
                %% check if solution has been found or how it changed if its not sufficient
                [ is_possible_solution , lambda2_smaller1_new] = check_IPS_solution(lambda_1, lambda_2, lambda_3, tolerance);
                
                % earlier break from while loop so that g does not get
                % changed any more
                if is_possible_solution
                    break
                end
                %            lambda2_old = lambda_2;
                
                % change the search direction and reduce step intervall if lambda2
                % passes one but is not in the required precision range.
                if lambda2_smaller1 ~= lambda2_smaller1_new
                    delta_g = - 0.5 * delta_g;              % Einbau intelligenter Schrittweitensteuerung wenn kein Fortschritt - haben es versucht, sind gescheitert... added break
                    %error('passed 1...')
                    lambda2_smaller1 = lambda2_smaller1_new;
                end
                             
%                 if abs( g - 17.253552526231005 ) < 1.e-15
%                     is_possible_solution;
%                     lambda2_smaller1;
%                     lambda2_smaller1_new
%                     x = 1;
%                 end
                                
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
%                 det_F = det(F)
%                 det_R = det(R)
%                 det_S = det(S)
%                 det_S_mirror = det(S_mirror)
%                 calculate invariant plane vector n_i etc.

                 % here the determinant of F changes for the second solution... 
                 % the first is the second of Khachaturyan...

                [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F,I,tolerance); %_kachaturyan2(F);

                % Note habit plane solutions come in pairs!
                
                isol = isol + 2; % increase counter for number of solutions found
                eps_s1 = slip_planes_between_burgerssteps( ds(is1,:), g, ns(is1,:), 'cubic');
                eps_s2 = slip_planes_between_burgerssteps( ds(is2,:), g, ns(is2,:), 'cubic');
                eps_s = [eps_s1, eps_s2];
                d = [ds_product(is1,:); ds_product(is2,:)];
                n = [ns_product(is1,:); ns_product(is2,:)];
                if mod(isol,100) == 0
                    isol
                end
                
                % Create Slip_solution objects and append them to object array 
                solutions.array( isol-1 ) =  Slip_solution(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, eps_s, d, n, m_aust' );
                solutions.array( isol )   =  Slip_solution(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, eps_s, d, n, m_aust' );                
            end
            
        end % end of loop for second slip system
%         if isol == 6
%             break
%        end
    end % end of loop for first slip system
%     if isol == 6
%         break
%    end
    
end % end of loop over considered mirror planes in martensite

fprintf('Total number of solutions for lambda_2 = 1 found is: n_sol = %i :\n', isol)

end


