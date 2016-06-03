% calculate initial eigenvalues without shear modification
% Ehl: ? function should return the three eigenvectors ? 1st eigenvector is
% was saved to lambda2_smaller1 --> what is it used for? --> changed, so that
% eigenvectors are saved now (to ev_1 ev_2, ev_3)
[ lambda_1, lambda_2, lambda_3, e1, e2, e3] = sorted_eig_vals_and_vecs( B3 );

epsilon = 1.e-9;
g_min = 9.5; % minimum number of atom layers before a step due to the (continuum) applied shear occurs (LIS)
g_initial = 300.0;
[isSolution, lambda2Bain_smallerl] = check_solution( lambda_1, lambda_2, lambda_3, epsilon);

for im =1:1 %size(m,1) % number of considered mirror planes in martensite
    
    m(im,:)
    
    m_mart = m(im,:); % mirror plane in martensite
    m_aust = inverse(cp) * m_mart';
    
    for is1 = 1:(size(ds,1)-1) % loop for first slip system
        for is2 = 2:(is1+1):size(ds,1) % loop for second one
            
            delta_g = 10.;
            g = g_initial;
            
            % first system
            d1 = cp * ds(is1,:)';
            n1 = inverse(cp) * ns(is1,:)';
            % Ehl: mirror_by_plane recieves the plane normal as 1st argument
            % and the vector to be mirrored as 2nd argument
            % --> therefore change the vectors in call of function?!
            d11 = mirror_by_plane(m_aust, d1, eye(3));
            n11 = mirror_by_plane(m_aust, n1, eye(3));
            % second system
            d2 = cp * ds(is2,:)';
            n2 = inverse(cp) * ns(is2,:)';
            % Ehl: mirror_by_plane recieves the plane normal as 1st argument
            % and the vector to be mirrored as 2nd argument
            % --> therefore change the vectors in call of function?!
            d22 = mirror_by_plane(m_aust, d2, eye(3));
            n22 = mirror_by_plane(m_aust, n2, eye(3));
            
            % S1 and S2 are shears related by mirror symmetry
            S1  = d1  * n1';
            S11 = d11 * n11';
            S2  = d2  * n2';
            S22 = d22 * n22';          
            
            % for loop over g (shear magnitude - m in Paper) values
            % write input function for range of g values (g_min, g_max) and
            % epsilon
            % e.g. g = 10.0 : 0.1 : 50.0  % maximum shear all 10 layers a step minmum all 50
            
            lambda2_smaller1 = lambda2Bain_smallerl;
            while ( ~isSolution && (g > g_min) )
                
                if g > g_initial
                    error('this should not happen...')
                end
                %m = sym('m');
                % symbolic inv() is possible
                
                % calculate first double shear matrix
                S = eye(3) + (1./g)*(S1 + S11);
                % calculate the sheared mirrorplane
                m_aust_sheared = inverse( S )* m_mart';
                
                % calculate the rotation of the mirror plane vector due
                % to the shear
                R = max_shear_rotation(m_aust, S);
                % Construct the net deformation gradient from the two
                % sheared sides
                F = 0.5*( inverse(R)*S + R*(eye(3) + (1./g)*(S2 + S22)) ) * B3; % composite block deformations
                
                % cf = sym( det( F - eye(3) ) ) == 0;
                % solutions = double( solve( cf, g) );   %, 'Real', true);
                
                % get new results
                % Ehl: in rank_one_kachaturyan wird F2 = F'*F verwendet ... so auch in der Veroeff.
                % rank_one_kachaturyan liefert dann andere Eigenwerte und läuft auf einen Fehler, weil der mittlere Eigenwert nicht nahe genug bei 1.0
                % deswegen hier F2 = F'*F ?
                F2 = (F'*F);
                [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F2 );
                
                %% check if solution has been found or how it changed if its not sufficient
                [ isSolution , lambda2_smaller1_new] = check_solution(lambda_1, lambda_2, lambda_3, epsilon);
                
                g;
                lambda_2;
                
                % change the search direction and reduce step intervall if lambda2
                % passes one but is not in the required precision range.
                % Kind of line search algorithm...
                if lambda2_smaller1 ~= lambda2_smaller1_new
                    delta_g = (-1./2.) * delta_g;
                    %error('passed 1...')
                end
                lambda2_smaller1 = lambda2_smaller1_new;
              
                % change g value.
                g = g - delta_g;    

            end % end while
            
            if isSolution
                %[ y1 ] = sorted_eig_vals_and_vecs( F ) % , y2, y3, e1, e2, e3]
                isSolution = false;
                % Ehl: in order to distinguish between the upper epsilon as criterion for the search-algorithm
                % and epsilon_0 as the magnitude of the shear (see paper Kachaturyan)
                % --> change of the return value 
                [eps_0, a1, a2, n1, n2, Q1, Q2] = rank_one_kachaturyan( F );
%                 [epsilon, a1, a2, n1, n2, Q1, Q2] = rank_one(F, eye(3) )

                % Ehl formatted output of the results
                fprintf('==================================================\n')
                fprintf('= results of calculation for a slip-system that  =\n')
                fprintf('= provides an invariant plane (after kachaturyan)=\n')
                fprintf('==================================================\n')
                fprintf('magnitude of shear: \n epsilon_0 = %1.4f \n', eps_0)
                fprintf('normal to habit plane (unit vector) - 1st solution: \n n_1 = \n')
                disp (n1)
                fprintf('normal to habit plane (unit vector) - 2nd solution: \n n_2 = \n')
                disp (n2)
                fprintf('shear direction of transformation (unit vector) - 1st solution: \n l_1 = \n')
                disp (a1)
                fprintf('shear direction of transformation (unit vector) - 2nd solution: \n l_2 = \n')
                disp (a2)
                fprintf('rotation matrix (for invariant planar match with parent) - 1st solution: \n R_I1 = \n')
                disp(Q1)
                fprintf('rotation matrix (for invariant planar match with parent) - 2nd solution: \n R_I2 = \n')
                disp(Q2)
                fprintf('==================================================\n')
                
                eps_0;
    
            end
%             %% if the solution for g is very high or low respectively, do not consider it
%             solution = solutions;
%             if solution > 10
%                 shear_amount(nr,1) = 7;
%             elseif solution < 10
%                 if solution < -10
%                     shear_amount(nr,1) = -7;
%                 elseif solution > - 10
%                     shear_amount(nr,1) = solution;
%                 end
%             end
%             shear_amount(nr,2) = nr;
            
            
        end
    end
    
end

