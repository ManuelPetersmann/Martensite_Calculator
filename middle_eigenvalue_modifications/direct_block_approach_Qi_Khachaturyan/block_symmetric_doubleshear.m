function [solutions] = block_symmetric_doubleshear(martensite, austenite)
% direct block approach of highly dislocated, blocky lath martesite
% microstructure after Qi,Khachaturyan,Morris 2014 Acta
% All calulations are carried out in the coordinate system of the parent phase
% returns object array of solutions for IPSs.

if (nargin < 2 && martensite.considered_plasticity == 1)
    austenite = 1; % set austenite to some random value to make function callable with just one argument
end

% specify type of solution array
martensite.IPS_solutions.array = Slip_solution();
% set calcuation method property in solution_array object
calculation_method = 'direct block approach, mirrorsym. & equal double-shears';
martensite.IPS_solutions.calculation_method = calculation_method;  
% create shorthand notation
solutions = martensite.IPS_solutions;

%% set numerical parameters und create solution object
numerical_parameters;

%% transform product phase slip systems to parent phase and combine all in one array
% assemble all shear directions, planes and dyads
[ds, ns, S, slip_combinations] = shear_dyads(martensite, austenite, true); % assemble miller- shear_dyads
solutions.slip_combinations = slip_combinations; % nr of possibilites nchoosek (k=2)

%% calculate only initial eigenvalues without shear modification to determine
% the direction from which side lambda2 = 1 is approached
[ lambda_1, lambda_2, lambda_3] = sorted_eig_vals_and_vecs( martensite.U'*martensite.U );
[~, lambda2_smaller1_initial] = check_IPS_solution( lambda_1, lambda_2, lambda_3, tolerance);

% lambda2_old = lambda_2;  % not used since it is not optimized incrementally...

isol = 0;
%% loop over mirror planes and slip systems
for im = 1:size(martensite.mirror_planes,1) % number of considered mirror planes in martensite
    
    m_mart = martensite.mirror_planes(im,:); % mirror plane in martensite
    m_aust = inverse(martensite.cp)' * m_mart'; % transformed plane in austenite
    
    % loop over slip system combinations
    for is1 = 1:(size(ds,1)-1) % loop for first slip system
        for is2 = (is1+1):size(ds,1) % loop for second one
            d1_mirr = mirror_vec_by_plane(m_aust, ds(is1,1:3), I);
            n1_mirr = mirror_vec_by_plane(m_aust, ns(is1,1:3), I);
            d2_mirr = mirror_vec_by_plane(m_aust, ds(is2,1:3), I);
            n2_mirr = mirror_vec_by_plane(m_aust, ds(is2,1:3), I);
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
                
                S1 =  I + (1./g)* (S(:,:,is1) + S(:,:,is2));
                S_mirror = I + (1./g)* (S11 + S22);
                
                % calculate the rotation of the mirror plane vector due to shear S1
                R = max_shear_rotation( m_aust, S1);
                
                % Construct the net deformation gradient from the two
                % sheared sides
                F = 0.5*( R*S1 + inverse(R)*S_mirror ) * martensite.U; % composite block deformations, matrix multiplication is distributive
                
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

                [y1,y3, d1, d2, h1, h2, Q1, Q2] = rank_one(F, I, tolerance );

                % Note habit plane solutions come in pairs!
                
                isol = isol + 2; % increase counter for number of solutions found
                eps_s1 = slip_planes_between_burgerssteps( ds(is1,1:3), g, ns(is1,1:3), 'cubic');
                eps_s2 = slip_planes_between_burgerssteps( ds(is2,1:3), g, ns(is2,1:3), 'cubic');
                eps_s = [eps_s1, eps_s2];
                d = [ds(is1,:); ds(is2,:)];
                n = [ns(is1,:); ns(is2,:)];
                if mod(isol,100) == 0
                    isol
                end
                % Create Slip_solution objects and append them to object array 
                solutions.array( isol-1 ) =  Slip_solution(F, I, y1, y3, d1, h1, Q1, Q1*martensite.U, eps_s, d, n, m_aust' );
                solutions.array( isol )   =  Slip_solution(F, I, y1, y3, d2, h2, Q2, Q2*martensite.U, eps_s, d, n, m_aust' );   
                solutions.array( isol-1 ).id = isol-1;
                solutions.array( isol ).id = isol;
            end
            
        end % end of loop for second slip system
    end % end of loop for first slip system

    
end % end of loop over considered mirror planes in martensite

if isol > 0 
disp(['number of potential solutions found = ', num2str(isol)])
solutions.solutions_available = 1;
end

end


