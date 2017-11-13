clc
clear all;

a_aust = 3.6017264; % for 140 Grad Celsius, 3.5975576 for 80 Grad Celsius
a_mart = 2.8807346; % for 140 Grad Celsius, 2.8790068 for 80 Grad Celsius- check if something changes 

Bain_and_Correspondence;

%% assemble slip systems in alpha
% since the shear is a substantial part of the transformation only 
% shear systems which are favorable in the b.c.c. lattices are considered. 
% the plane and direction families are {110}_alpha, {112}_alpha,
% <111>_alpha, <110>_alpha
plane_families_bcc =     [ [1 1 0]
                           [1 1 2] ];   % must be written with linebreak or ";" between vectors!                     
direction_families_bcc = [ [1 1 1]
                           [1 1 0] ];
                       
count_directions_extra = true;
                       
% find all possible combination (including different shear directions)
[martensite.slip_planes, martensite.slip_directions] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );
%[ ns_product, ds_product ] = independent_slipsystems( plane_families_bcc, direction_families_bcc, count_directions_extra );

plane_families_fcc =     [ [1 1 1] ];
direction_families_fcc = [ [1 1 0]; [1 1 2] ];
[austenite.slip_planes, austenite.slip_directions] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);
%[ ns_parent, ds_parent] = independent_slipsystems(plane_families_fcc,direction_families_fcc,count_directions_extra);

martensite.considered_plasticity = 3; % 1-mart, 2-aust, 3-both mart and aust slip systems
    
cpps_gamma = all_from_family_perms( [1 1 1] );
austenite.CPPs = cpps_gamma;

[ds, ns, S, slip_combinations] = shear_dyads(martensite, austenite, false); % assemble normed- shear_dyads

% 'Kurdjumov Sachs directions [110]_aust || [111]_mart';
% densest packed direction in austenite
% KS = u !!!!
us = all_from_family_perms( [1 1 0] ); %, false ); % second argument sorts out sign-ambiguous vectors, i.e. [1 1 0] = [-1 -1 0]
us = us / sqrt(2);


u = us(1,:)';
%B3 * u
%%

numerical_parameters;

%length(S)
%u = [1. 1. 1.]'

% fid = fopen('equivalent_screws','w');
% fprintf(fid,'%s \t \t %s \t \t %s',' s ',' b ',' m');

counter = 1;
for iu = 1:length(u)
    u = us(iu,:)';
    u2  = B3 * u;
    R   = rotation_between_vectors( u2,     u );
    %try
    [angle, ax] = rotmat_to_axis_angle( R );
    angle_mod = angle;
    ax_mod = ax;
    %catch
    % it has been tested that the specific directions of <110> that
    % remain unrotated by B3 cannot be made any better using shears!
    %end
    u2_unrot = R * u2;
    %
    for is1 = 1:length(S)
        for is2 = is1+1:length(S)
            
            %        [screw_syst] = equivalent_shear_screwdisloc( ds(i,1:3), ns(i,1:3) );
            %        fprintf(fid,'\n %s \t \t %s \t \t %s', mat2str(ds(i,1:3)), mat2str(screw_syst(1:3)), mat2str(screw_syst(4:6)) );
            
            res_old = norm( u - u2_unrot ); % If an invariant line is specified a -
            % priori and it does not ly on the cone of undistorted line
            % then, even it is rotated back to its inital position it may be
            % distorted. Also the rotation to bring it back may be quite large. In
            % the following it is tried to minimize the rotation and make a given
            % line fully undistorted by adding additional plastic shears to the
            % total deformation of individual martensite laths
            % res... residual - if it is zero then the vector is fully invariant
            
            delta_eps = delta_eps_initial;
            eps1 = eps_initial;
            eps2 = eps_initial;
            % g = g_initial / 1./(norm(d11)*norm(n11));
            is_possible_solution = false;
            S_accummulated = eye(3);
            
            while ( ( res_old > tolerance )  ...
                 && (eps1 < eps_max)        && (eps2 < eps_max) )  
                % if the solution for 'eps' is very high or low respectively, do not consider it
                if ((eps1 < eps_initial) || (eps2 < eps_initial))
                    error('this should not happen - fix code...')
                end
                
                u2_dS1 = B3 *  S_accummulated * ( eye(3) + delta_eps*S(:,:,is1) ) * u;
                %try
                % calculate R such that u2_dS1 is unrotated
                R_dS1 = rotation_between_vectors( u2_dS1, u );
                %catch
                %    continue
                %end
                u2_dS1_unrot = R_dS1 * u2_dS1;
                new_res_dS1 = norm( u - u2_dS1_unrot ); % check deviation from invariance
                
                
                u2_dS2 = B3 *  S_accummulated * ( eye(3) + delta_eps*S(:,:,is2) ) * u;
                %try
                % calculate R such that u2_dS1 is unrotated
                R_dS2 = rotation_between_vectors( u2_dS2, u );
                %catch
                %    continue
                %end
                u2_dS2_unrot = R_dS2 * u2_dS2;
                new_res_dS2 = norm( u - u2_dS2_unrot ); % check deviation from invariance
                
                
                % choose the shear that approached a solution quicker with the
                % SAME shear magnitude  - delta_eps
                if(new_res_dS1 < new_res_dS2)
                    % here generally the minimum should be taken if the solution has not been passed already
                    % from one search direction, otherwise the shear amplitude is adopted first and only then it is
                    % checked again wheter a new minimum does not overshoot the solution from one direction
                    if ( res_old < new_res_dS1 ) % ( ( lambda2_smaller1_shear1 ~= lambda2_smaller1_initial ) ||
                        % i.e. the solution is passed or delta_lambda2_to_1 does not decrease, which should not happen - cut back
                        delta_eps = stepwidth_change * delta_eps;
                    else
                        res_old = new_res_dS1;
                        eps1 = eps1 + delta_eps;
                        S_accummulated = S_accummulated * ( I + delta_eps * S(:,:,is1) );
                        R_mod = R_dS1;
                    end
                else % new_res_dS1 > new_res_dS2)
                    if ( res_old < new_res_dS2 ) % ( ( lambda2_smaller1_shear2 ~= lambda2_smaller1_initial ) ||
                        delta_eps = stepwidth_change * delta_eps;
                    else
                        res_old = new_res_dS2;
                        eps2 = eps2 + delta_eps;
                        S_accummulated = S_accummulated * ( I + delta_eps * S(:,:,is2) );
                        R_mod = R_dS2;
                    end
                end
                
                if (delta_eps < tolerance)
                    break
                end
                
            end % end while
            
            
            if ( res_old < tolerance )
                F_tot = R_mod * B3 *  S_accummulated;
                [U_toal , R_total] = polardecomposition( F_tot );
                % this is executed for constant inv line u  before: [angle, ax] = rotmat_to_axis_angle( R );
                [angle_mod, ax_mod] = rotmat_to_axis_angle( R_total );
                
                if (angle_mod < 6.)
                    [angle_R, ax_R] = rotmat_to_axis_angle( R_mod );
                    
                    if is_rank_one_connected(F_tot,eye(3), 1.e-2)
                        
                        ILS_laths(:,:,counter) = F_tot;
%                         is1
%                         is2
%                         res = norm( u - u2_unrot )
%                         res_old % modified res
%                         ax
%                         ax_R
%                         ax_mod
%                         angle
%                         angle_R
%                         angle_mod
%                         eps1
%                         eps2
%                         
%                         isol = isol + 2; % increase counter for number of solutions found
%                         if mod(isol,500)==0
%                             isol
%                             %pause(1);
%                         end
%                         eps_s = [eps1; eps2];
%                         d = [ds(is1,:); ds(is2,:)];
%                         n = [ns(is1,:); ns(is2,:)];
%                         
%                         % Create Slip_solution objects and append them to object array;
%                         % PET 10.10.17: replaced 'isol' and 'eps' wit y1 and y2
%                         solutions.array( isol-1 ) =  Slip_solution(F, I, y1, y3, d1, h1, Q1, Q1*martensite.U, eps_s, d, n );
%                         solutions.array( isol )   =  Slip_solution(F, I, y1, y3, d2, h2, Q2, Q2*martensite.U, eps_s, d ,n );
%                         % reduced contructor could look like
%                         %  solutions.array( isol-1 ) =  Slip_solution(F, I, martensite.U, eps_s, d, n );
%                         solutions.array( isol-1 ).id = isol-1;
%                         solutions.array( isol ).id = isol;
                        counter = counter + 1;
                    end
                end
            end
            
            
            
        end
    end
end
counter = counter-1;
counter
