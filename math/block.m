
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
                
                % for loop over m... values
                
                % S1 and S2 are shears related by mirror symmetry
                m = sym('m');
                S1  = (1./m)*(d1  * n1');
                S11 = (1./m)*(d11 * n11');
                S2  = (1./m)*(d2  * n2');
                S22 = (1./m)*(d22 * n22');
                % S1 = get_shear_d_n( d1, n1);
                % S2 = get_shear_d_n( d2, n2);
                
                S = eye(3) + S1 + S11;
                m_aust_sheared = inverse( S )* m_mart'; % symbolic inv() is possible
                
                R = max_shear_rotation(m_aust, S);
                
                F = sym( 0.5*( inverse(R)*S + R*(eye(3) + S2 + S22) * B3 ) ); % composite block deformations
                
                cf = sym( det( F - eye(3) ) ) == 0;
                
                solutions = double( solve( cf, g) );   %, 'Real', true);
                
                % If for this equation no solution could be found assign 8
                if isempty( solutions )
                    shear_amount( nr,1) = 8;
                    shear_amount( nr,2) = i;
                    nr = nr +1;
                    continue
                end
                % if the solution for g is very high or low respectively set it to 7
                % or -7 for a good readability of the results
                %     for si = 1:2
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
                
                %end
        end
    end
    
end

