function [solutions] = block_symmetric_doubleshear(B, cp, ms, ns, ds )
% call: block_symmetric_doubleshear(B, cp, ms, ns, ds)
% B... Bain strain, cp - B*Correspondance matrix, ms...mirror planes in
% alpha', ns...slip system normals in alpha, ds... slip directios in alpha
% returns object array of solutions for IPSs.

%% initalize some other vars
solutions = Solution_array( Slip_solution() ); % Construct array with type of solution -> After this line, Solution_array.array is no longer a double 
epsilon = 1.e-11; % accuracy for middle valued eigenvalue
g_min = 4.;
g_initial = 100.0; 
delta_g_initial = 5.; % must be set reasonably so that a solution is found i.e. g_min not to high, delta_g not to high either
% here it reaches as lowest g=5 with delta_g = 5. everything lower will be disregarded
I = eye(3);
isol = 0; % counter for number of solutions

%k = 0.5 % stepwidth scaling factor



%% calculate only initial eigenvalues without shear modification to determine
% the direction from which side lambda2 = 1 is approached
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( B'*B );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, epsilon);

lambda2_old = lambda_2;

%% loop over mirror planes and slip systems
for im = 1:size(ms,1) % number of considered mirror planes in martensite
    
    m_mart = ms(im,:); % mirror plane in martensite
    m_aust = inverse(cp)' * m_mart'; % transformed plane in austenite
    
    % loop over slip system combinations
    for is1 = 1:(size(ds,1)-1) % loop for first slip system
        for is2 = (is1+1):size(ds,1) % loop for second one
            
            % ############  transform bcc systems to austenite ############
            % first system
            d1 = cp * ds(is1,:)'; % martensite.vec_from_coords( ds(is1,:) )'; - not necessary for cubic lattice
            n1 = inverse(cp)' * ns(is1,:)';
            d11 = mirror_by_plane(m_aust, d1, I);
            n11 = mirror_by_plane(m_aust, n1, I);
            % second system
            d2 = cp * ds(is2,:)';
            n2 = inverse(cp)' * ns(is2,:)';
            d22 = mirror_by_plane(m_aust, d2, I);
            n22 = mirror_by_plane(m_aust, n2, I);
            % S_i and S_ii are shears related by mirror symmetry
            
            % display('normed run');
%             S1  = (d1  * n1') * 1./(norm(d1)*norm(n1));
%             S2  = (d2  * n2') * 1./(norm(d2)*norm(n2));
%             S11 = (d11 * n11') * 1./(norm(d11)*norm(n11));
%             S22 = (d22 * n22') * 1./(norm(d22)*norm(n22));
            
            % display('Khachaturyan run');
            S1  = (d1  * n1') ;
            S2  = (d2  * n2') ;
            S11 = (d11 * n11') ;
            S22 = (d22 * n22') ;
            
            %% modify shear value in Blocks until lambda2 = 1
            delta_g = delta_g_initial;
            g = g_initial;
            % g = g_initial / 1./(norm(d11)*norm(n11));
            is_possible_solution = false;
            lambda2_smaller1 = lambda2_smaller1_initial;
            while ( ~ is_possible_solution && (g > g_min) )  % if the solution for g is very high or low respectively, do not consider it
                
                if g > g_initial
                    error('this should not happen - fix code...')
                end
                
                % given the slip direction d and the slip plane normal n and the shear magnitude
                % 1/g (i.e. the lattice has a step after each g planes after the shear)
                % calculate double shear matrizes
                
                % shears are not commutative: Sx = (I + S1)*(I + S2) \uneq Sy = (I + S2)*(I + S1)
                % S = (I + (1./g)*S1)*(I + (1./g)*S2)
                % S_mirror = (I + (1./g)*S11)*(I + (1./g)*S22)
                % verified that like Khachaturyan writes it, it is the same
                % as the above multiplied version, if the first shear comes first
                S =  I + (1./g)* (S1 + S2);
                S_mirror = I + (1./g)* (S11 + S22);
                
                % calculate the rotation of the mirror plane vector due to the shear
                R = max_shear_rotation( m_aust, S);
                
                % Construct the net deformation gradient from the two
                % sheared sides
                F = 0.5*( R*S + inverse(R)*S_mirror ) * B; % composite block deformations, matrix multiplication is distributive
                
                % get new results
                [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
                
                %% check if solution has been found or how it changed if its not sufficient
                [ is_possible_solution , lambda2_smaller1_new] = check_IPS_solution(lambda_1, lambda_2, lambda_3, epsilon);
                
                if(abs(lambda2_old - lambda_2) < 1.e-15)
                    break
                end
                
                % change the search direction and reduce step intervall if lambda2
                % passes one but is not in the required precision range.
                if lambda2_smaller1 ~= lambda2_smaller1_new
                    delta_g = - 0.5 * delta_g;              % Einbau intelligenter Schrittweitensteuerung wenn kein Fortschritt - haben es versucht, sind gescheitert... added break
                    %error('passed 1...')
                end
                
                
%                 if abs( g - 17.253552526231005 ) < 1.e-15
%                     is_possible_solution;
%                     lambda2_smaller1;
%                     lambda2_smaller1_new
%                     x = 1;
%                 end
                
                lambda2_smaller1 = lambda2_smaller1_new;
                
                lambda2_old = lambda_2;
                
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
%                [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F, I ) % here the determinant of F changes for the second solution... the first is the second of Khachaturyan...
%                 Q1
%                 a1
%                 h1
%                 Q2
%                 a2
%                 h2
%                 det_ST1 = det( Q1*F)
%                 det_ST2 = det( Q2*F)
                [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F,I);
%                 Q1
%                 a1
%                 h1
%                 Q2
%                 a2
%                 h2
                
%                 det_LT1 = det( Q1*B)
%                 det_LT2 = det( Q2*B)
%                 det_ST1 = det( Q1*F)
%                 det_ST2 = det( Q2*F)
%                 a1
%                 a22
%                 a2
%                 a11
%                 h1
%                 h22
%                 h2
%                 h11

                % Note habit plane solutions come in pairs!
                
                isol = isol + 2; % increase counter for number of solutions found
                
                % Create Slip_solution objects and append them to object array 
                solutions.array( isol-1 ) =  Slip_solution(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, g, ds(is1,:), ns(is1,:), ds(is2,:), ns(is2,:), m_aust' );
                solutions.array( isol )   =  Slip_solution(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, g, ds(is1,:), ns(is1,:), ds(is2,:), ns(is2,:), m_aust' );
                

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


fprintf('number of potential solutions found: n_sol = %i :\n', isol)

end


