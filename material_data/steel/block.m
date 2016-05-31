% calculate initial eigenvalues without shear modification
% Ehl: ? function should return the three eigenvectors ? 1st eigenvector is
% now saved to lambda2_smaller1 --> what is it used for? --> changed, so that
% eigenvectors are saved now (to ev_1 ev_2, ev_3)
[ lambda_1, lambda_2, lambda_3, ev_1, ev_2, ev_3] = sorted_eig_vals_and_vecs( B3 );

% epsilon is a criterion to check, if the transformation provides an
% invariant plane --> eigenvalue lambda_2 ~= 1.0
epsilon = 1.e-6;
% Ehl: what is this value for?
g_max = 50.0;

[noSolution, lambda2_smaller1] = check_solution( lambda_1, lambda_2, lambda_3, epsilon);

for im = 1:size(m,1) % number of considered mirror planes in martensite
    
    m(im,:)
    
    m_mart = m(im,:); % mirror plane in martensite
    % express mirror plane from bcc-lattice(martensite) in fcc-lattice(austenite)
    % see e.g. Qi. et al (2014) Eq.(22)
    m_aust = inverse(cp) * m_mart'; 
    
    for is1 = 1:1 %(size(ds,1)-1) % loop for first slip system
        for is2 = 2:2 %(is1+1):size(ds,1) % loop for second one
            
            delta_g = 20.;
            g = 5.0;
            
            % first system
            d1 = cp * ds(is1,:)';
            n1 = inverse(cp) * ns(is1,:)';
            % Ehl: ?changed orde of d1/n1 and m_aust, in order to fit to the definition of the function mirror_by_plane 
            % --> normal of the mirror plane comes in first ?
            d11 = mirror_by_plane(m_aust, d1, eye(3));
            n11 = mirror_by_plane(m_aust, n1, eye(3));
            % second system
            d2 = cp * ds(is2,:)';
            n2 = inverse(cp) * ns(is2,:)';
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
            
            while ( noSolution && (g < g_max) )
                %m = sym('m');
                % symbolic inv() is possible
                
                % calculated first double shear matrix
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
                [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F );
                
                %% check if solution has been found or how it changed if its not sufficient
                [ isSolution , lambda2_smaller1_new] = check_solution(lambda_1, lambda_2, lambda_3, epsilon);
                
                % change the search direction and reduce step intervall if lambda2
                % passes one but is not in the required precision range.
                % Kind of line search algorithm...
                if lambda2_smaller1 ~= lambda2_smaller1_new
                    delta_g = (-1./2.) * delta_g;
                    display('lololololololol')
                end
                lambda2_smaller1_new = lambda2_smaller1;
              
                g
                lambda_2
                
                % change g value.
                g = g + delta_g;    

            end % end while
            
            if noSolution == false
                sorted_eig_vals_and_vecs( F )
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

