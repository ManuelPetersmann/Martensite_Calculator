
for im =1:size(m,1) % number of considered mirror planes in martensite
    
    m_mart = m(im,:); % mirror plane in martensite
    m_aust = inverse(cp) * m_mart';
    
    for is1 = 1:(size(ds,1)-1) % loop for first slip system
        for is2 = (is1+1):size(ds,1) % loop for second one
            %for m = 10.0 : 0.1 : 30.0  % maximum shear all 10 layers minmum all 30 
                
                d1 = cp * ds(is1,:)';
                n1 = inverse(cp) * ns(is1,:)';
                d11 = mirror_by_plane(d1, m_aust, eye(3));
                n11 = mirror_by_plane(n1, m_aust, eye(3));
                % second system
                d2 = cp * ds(is2,:)';
                n2 = inverse(cp) * ns(is2,:)';
                d22 = mirror_by_plane(d1, m_aust, eye(3));
                n22 = mirror_by_plane(n1, m_aust, eye(3));
                
                % for loop over g (shear magnitude - m in Paper) values
                % write input function for range of g values (g_min, g_max) and
                % epsilon
                
                %% S1 and S2 are shears related by mirror symmetry
                while ((lambda_2 - 1.) > 1.e-6) || (g > g_max)  %epsilon
                    %m = sym('m');
                    S1  = (1./g)*(d1  * n1');
                    S11 = (1./g)*(d11 * n11');
                    S2  = (1./g)*(d2  * n2');
                    S22 = (1./g)*(d22 * n22');
                    % S1 = get_shear_d_n( d1, n1);
                    % S2 = get_shear_d_n( d2, n2);
                    
                    S = eye(3) + S1 + S11;
                    m_aust_sheared = inverse( S )* m_mart'; % symbolic inv() is possible
                    
                    R = max_shear_rotation(m_aust, S);
                    
                    F = sym( 0.5*( inverse(R)*S + R*(eye(3) + S2 + S22) * B3 ) ); % composite block deformations
                    
                    cf = sym( det( F - eye(3) ) ) == 0;
                    
                    solutions = double( solve( cf, g) );   %, 'Real', true);
                    
                    %% if there is a solution check if the other eigenvalues
                    % straddle 1! I forgot that last time...
                    
                    % if the solution for g is very high or low respectively, do not consider it
                    solution = solutions;
                    if solution > 10
                        shear_amount(nr,1) = 7;
                    elseif solution < 10
                        if solution < -10
                            shear_amount(nr,1) = -7;
                        elseif solution > - 10
                            shear_amount(nr,1) = solution;
                        end
                    end
                    shear_amount(nr,2) = nr;
                    
                    %% if there is no solution
                    if isempty( solutions )
                        g = g + delta_g; 
                        continue
                    end
                    % increase m value. Make kind of line search algo to
                    % find how much it can increase so that not everything
                    % is calculated...
                    
                end % end while
        end
    end
    
end

