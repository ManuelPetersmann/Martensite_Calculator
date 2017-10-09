function [solutions] = interface_defects_doubleshear_MarescaCurtin(B, cp, ns_P2, ds_P2, ns_P3, ds_P3, phi)
% call: interface_defects_doubleshear_MarescaCurtin(B, cp, ms, ns, ds )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of shape deformation following:
% F. Maresca, W.A. Curtin, The austenite/lath martensite interface in steels:
% Structure, athermal motion, and in-situ transformation strain revealed by simulation and theory, 
% Acta Materialia (2017), doi: 10.1016/j.actamat.2017.05.044.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables:
% B    - Bain strain, 
% cp   - B*Correspondance matrix, 
% ns_2 - slip planes for 1st lattice defect --> P^(2)f
% ds_2 - slip directions for 1st lattice defect --> P^(2)
% ns_3 - slip planes for 2nd lattice defect --> P^(3)
% ds_3 - slip directions for 2nd lattice defect --> P^(3)
% phi  - prescribed orientation relationship
%
%
% Output variables:
% solutions - object array of solutions for IPSs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solutions = Solution_array( Slip_solution ); % Construct array with type of solution -> After this line, Solution_array.array is no longer a double 

numerical_parameters;
beta_min = 1.; % average step heigth of a physical interface must be greater than one fcc planar spacing
beta_initial = 20.; 
delta_beta_initial = 0.25; % must be set reasonably so that a solution is found
R_psi = zeros(3); % initialization
R_phi = zeros(3); % initialization

% precalculations for P^(2) - done here as m2 is constant for prescribed OR phi
zeta = acosd(1/3); % angle between (-1 0 1)_bcc and (011)_bcc
m2 = (2*sind(phi))/(sqrt(3)*sind(zeta)*sind(zeta-phi)); % magnitude of shear P^(2)

% calculate homogeneous fcc to bcc lattice transformation S
% psi = acosd(dot([0 0 1]_fcc,[1 1 1]_fcc)/(norm([0 0 1]_fcc)*norm([1 1 1]_fcc))) - 45°
%     = 54.7356°-45° = 9.7356°
% see also Maresca & Curtin (2017) Fig.2
psi = 9.7356; % angle psi for rotation around [100]_bcc which aligns (111)_fcc||(011)_bcc

% omega = phi_NW-phi
omega = 5.26-phi;

R_psi = vrrotvec2mat([1  0 0 (psi*pi/180.)]); % rotation R_psi by the angle psi around [100]_bcc aligns (111)_fcc||(011)_bcc
R_phi = vrrotvec2mat([1  1 1 (omega*pi/180.)]); % rotation R_phi by angle omega = phi_NW - phi around [011]_bcc||[111]_fcc

S = R_phi * R_psi * B; 


%% calculate only initial eigenvalues without shear modification to determine
% the direction from which side lambda2 = 1 is approached
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( B'*B );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, tolerance);

lambda2_old = lambda_2;
    
% loop over slip system combinations: 
% each case one slip from P^(2) + one from P^(3)
for is1 = 1:size(ds_P2,1) % loop over all slip systems for P^(2)
    for is2 = 1:size(ds_P3,1) % loop over all slip systems for P^(3)
            
        % ############  transform bcc systems to fcc ############
        % first system
        d_P2 = cp * ds_P2(is1,:)'; % martensite.vec_from_coords( ds(is1,:) )'; - not necessary for cubic lattice
        n_P2 = inverse(cp)' * ns_P2(is1,:)';

        % second system - already in fcc, don't have to be transformed
        d_P3 = ds_P3(is2,:)';
        n_P3 = ns_P3(is2,:)';
        % d_P3 = cp * ds_P3(is2,:)';
        % n_P3 = inverse(cp)' * ns_P3(is2,:)';
        
        % display('normed run');
        %             S1  = (d1  * n1') * 1./(norm(d1)*norm(n1));
        %             S2  = (d2  * n2') * 1./(norm(d2)*norm(n2));
                 
        %% modify shear value in Blocks until lambda2 = 1
        delta_beta = delta_beta_initial;
        beta = beta_initial;
        
        % g = g_initial / 1./(norm(d11)*norm(n11));
        
        is_possible_solution = false;
        lambda2_smaller1 = lambda2_smaller1_initial;
        while ( ~ is_possible_solution && (beta > beta_min) )  % if the solution for beta is very high or < 1 , do not consider it
            
            if beta > beta_initial
%                 error('this should not happen - fix code...')
                break
            end
            
            % given the slip directions d_P2, d_P3,
            % the slip plane normals n_P2, n_P3
            % the orientation relationship phi
            % and the average step height beta
            % calculate the shears P^(2) and P^(3) and with these the observed shear P^(1)
            
            % calculate P^(2) --- constant for given OR
            P2 = I + m2 * (d_P2*n_P2');
            
            % calculate P^(3)
            m3 = (1/beta)*sqrt(3/2); % magnitude of shear P^(3)
            P3 = I + m3 * (d_P3*n_P3');
            
            % Construct the net deformation gradient 
            % with help of P^(2), P^(3) and S ( P^(1) := F )
            %
            % F = R_{Delta} * S * P3 * P2
            %
            F = S * P3 * P2;
            
            % get new results
            [ lambda_1, lambda_2, lambda_3 ] = sorted_eig_vals_and_vecs( F'*F );
            
            %% check if solution has been found or how it changed if its not sufficient
            [ is_possible_solution , lambda2_smaller1_new] = check_IPS_solution(lambda_1, lambda_2, lambda_3, tolerance);
            
            if(abs(lambda2_old - lambda_2) < 1.e-15)
                break
            end
            
            % change the search direction and reduce step intervall if lambda2
            % passes one but is not in the required precision range.
            if lambda2_smaller1 ~= lambda2_smaller1_new
                delta_beta = - 0.5 * delta_beta;              % Einbau intelligenter Schrittweitensteuerung wenn kein Fortschritt - haben es versucht, sind gescheitert... added break
                %error('passed 1...')
            end
            
            lambda2_smaller1 = lambda2_smaller1_new;
            
            lambda2_old = lambda_2;
            
            % change value of beta.
            beta = beta - delta_beta;
            % find beta (average step height normalized by (a_fcc / sqrt(3))- in Paper Maresca & Curtin 2017)
            % within specified limits (beta_min=1,  beta_max)
            
            
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
            [eps_0, a1, a2, h1, h2, Q1, Q2] = rank_one(F,I,tolerance);
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
            solutions.array( isol-1 ) =  Slip_solution(F, I, isol-1, eps_0, a1, h1, Q1, Q1*B, 1./m2, ds_P2(is1,:), ns_P2(is1,:), beta/sqrt(1.5), ds_P3(is2,:), ns_P3(is2,:));
            solutions.array( isol )   =  Slip_solution(F, I, isol,   eps_0, a2, h2, Q2, Q2*B, 1./m2, ds_P2(is1,:), ns_P2(is1,:), beta/sqrt(1.5), ds_P3(is2,:), ns_P3(is2,:));
            
            
        end
        
    end % end of loop for second slip system
    %         if isol == 6
    %             break
    %        end
end % end of loop for first slip system
%     if isol == 6
%         break
%    end


fprintf('Total number of solutions for lambda_2 = 1 found is: n_sol = %i :\n', isol)

end


