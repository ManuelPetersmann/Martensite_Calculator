function solutions = invariant_line_strain(martensite, austenite)
% calculates invariant line strain (rototated, unstreched) - per default
% the invariant line is assummed to be the close packed direction in
% austenite!

if isempty( martensite.invariant_lines )
    us = austenite.CP_dirs;
    % VECTORS MUST BE NORMED!
end
us

%% NOTE: martensite is a handle class so everything that is set here is set everywhere!

% create shorthand notation
solutions = martensite.ILS_solutions;
% specify type of solution array
solutions.array = ILS_solution();
% set calcuation method property in solution_array object
solutions.calculation_method = 'variable doubleshear incremental optimization to give invariant line';


%% set numerical parameters (see file numerical_parameters.m)
numerical_parameters;

%% transform product phase slip systems to parent phase and combine all in one array
% assemble all shear dyads in austenite, array of directions, planes and in
% respective phase (miller indizes)
[ds, ns, S, slip_combinations] = shear_dyads(martensite, austenite, false); % assemble normed- shear_dyads
solutions.slip_combinations = slip_combinations; % nr of possibilites nchoosek (k=2)

% fid = fopen('equivalent_screws','w');
% fprintf(fid,'%s \t \t %s \t \t %s',' s ',' b ',' m');

isol = 0;
for iu = 1:length(us)
    u = us(iu,:)';
    u = u / norm(u); % VECTOR MUST BE NORMED!
    u2  = martensite.U * u;
    R_Bain   = rotation_between_vectors( u2,     u );
    %try
    %[angle_Bain, ax_Bain] = rotmat_to_axis_angle( R_Bain );
    %angle_inclusion = angle_Bain;
    %ax_inclusion = ax_Bain;
    %catch
    % it has been tested that the specific directions of <110> that
    % remain unrotated by B3 = martensite.U cannot be made any better using shears!
    %end
    u2_unrot = R_Bain * u2;
    %
    for is1 = 1:length(S)
        for is2 = is1+1:length(S)
            
            %  u = cross( ns(is1,1:3), ns(is2,1:3) );
            %  h = cross( ds(is1,1:3), ds(is2,1:3) );
            
            %  [screw_syst] = equivalent_shear_screwdisloc( ds(i,1:3), ns(i,1:3) );
            %  fprintf(fid,'\n %s \t \t %s \t \t %s', mat2str(ds(i,1:3)), mat2str(screw_syst(1:3)), mat2str(screw_syst(4:6)) );
            
            res_initial = norm( u - u2_unrot ); % If an invariant line is specified a -
            res_old = res_initial;
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
            S_accummulated = eye(3);
            %
            [new_res_dS1, ~] = perp_ILS( martensite.U, ( eye(3) + delta_eps*S(:,:,is1) ), u);
            [new_res_dS2, ~] = perp_ILS( martensite.U, ( eye(3) + delta_eps*S(:,:,is1) ), u);
            % stop immediately if the solution does not improve at all
            % this could be moved outside the while but then the code
            % above must be written twice...
            if ( res_old > new_res_dS1 ) && ( res_old > new_res_dS2 )
                
                while ( ( res_old > vec_residual )  ...
                        && (eps1 < eps_max)        && (eps2 < eps_max) )
                    % if the solution for 'eps' is very high or low respectively, do not consider it
                    if ((eps1 < eps_initial) || (eps2 < eps_initial))
                        error('this should not happen - fix code...')
                    end
                    
                    [new_res_dS1, R_dS1] = perp_ILS( martensite.U*S_accummulated, ( eye(3) + delta_eps*S(:,:,is1) ), u);
                    [new_res_dS2, R_dS2] = perp_ILS( martensite.U*S_accummulated, ( eye(3) + delta_eps*S(:,:,is1) ), u);
                    
                    %                 u2_dS1 = martensite.U *  S_accummulated * ( eye(3) + delta_eps*S(:,:,is1) ) * u;
                    %                 %try
                    %                 % calculate R such that u2_dS1 is unrotated
                    %                 R_dS1 = rotation_between_vectors( u2_dS1, u );
                    %                 %catch
                    %                 %    continue
                    %                 %end
                    %                 u2_dS1_unrot = R_dS1 * u2_dS1;
                    %                 new_res_dS1 = norm( u - u2_dS1_unrot ); % check deviation from invariance
                    %
                    %
                    %                 u2_dS2 = martensite.U *  S_accummulated * ( eye(3) + delta_eps*S(:,:,is2) ) * u;
                    %                 %try
                    %                 % calculate R such that u2_dS1 is unrotated
                    %                 R_dS2 = rotation_between_vectors( u2_dS2, u );
                    %                 %catch
                    %                 %    continue
                    %                 %end
                    %                 u2_dS2_unrot = R_dS2 * u2_dS2;
                    %                 new_res_dS2 = norm( u - u2_dS2_unrot ); % check deviation from invariance
                    
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
                    
                    if (delta_eps < delta_eps_tolerance) % 0.5^16 = 1.e-5||
                        break
                    end
                    
                end % end while
            end
            
            
            if ( res_old < vec_residual )
                
                F_tot = R_mod * martensite.U *  S_accummulated;
                LT    = R_mod * martensite.U;     
                
                isol = isol + 1;
                if mod(isol,100)==0
                    isol
                    %pause(1);
                end
                
                eps_s = [eps1; eps2];
                d = [ds(is1,:); ds(is2,:)];
                n = [ns(is1,:); ns(is2,:)];
                
                %solutions.array( isol-1 ) =  ILS_solution(u, F_tot, LT); % F_tot = ST...shape transformation
                %solutions.array( isol-1 ).slip = Slip_systems( eps_s, d, n );
                solutions.array( isol )      =  ILS_solution(u, F_tot, LT, R_Bain);
                solutions.array( isol ).slip = Slip_systems( eps_s, d, n );    
            %else
            %    neg_no_convergence_to_ILS = neg_no_convergence_to_ILS +1;
            end
                
        end
    end
end

%% SELECTION IS NOW DONE INDIVIDUALLY

%% selection criteria
% max_rot_angle_inclusion = 15. % degree
% %added_mass_angle_tolerance = rot_angle_tolerance;
% theta_CP_max = 2. % degree
% lambda2_ips_tolerance_lath = 2.e-2
% 
% neg_no_convergence_to_ILS = 0;
% neg_theta_CP = 0;
% neg_far_from_IPS = 0;
% neg_rot = 0;


%             [~ , R_total] = polardecomposition( F_tot );
%             [angle_inclusion, ax_inclusion] = rotmat_to_axis_angle( R_total );
%             added_mass_angle = angle_ax_u + angle_inclusion; % added_mass_angle XD
%             if (angle_inclusion < max_rot_angle_inclusion)
%                 %if (added_mass_angle < added_mass_angle_tolerance)
%                 %[angle_lattice, ax_lattice] = rotmat_to_axis_angle( R_mod );
%                 
%                 [bool, lambda2] = is_rank_one_connected(F_tot,eye(3), lambda2_ips_tolerance_lath);
%                 if bool
%                     
%                     LT = R_mod*martensite.U;
%                     
%                     [ theta_CP, closest_111aust_to_CP ] = min_misorientation( austenite.CPPs, LT, true );
%                     %theta_CP
%                     if theta_CP < theta_CP_max
%                         
%                         %                         res = norm( u - u2_unrot )
%                         %                         res_old % modified res
%                         
%                         %                         ax_Bain
%                         %                         ax_lattice
%                         %                         ax_inclusion
%                         
%                         %                         angle_Bain
%                         %                         angle_lattice
%                         %                         angle_inclusion
%                     else
%                         neg_theta_CP = neg_theta_CP +1;
%                     end
%                 else
%                     neg_far_from_IPS = neg_far_from_IPS +1;
%                 end
%             else
%                 neg_rot = neg_rot + 1;
%             end

% disp( [num2str(nchoosek(length(S),2)*size(us,1)),' combinations {u_i,S1_j,S2_k} tested, ', ...
%     num2str(isol),' solutions found, ' num2str(neg_no_convergence_to_ILS), ' did not converge to ILS'] );
% disp( ['First crit: rotation target function of lath > ', num2str(max_rot_angle_inclusion),' degree --> ', num2str(neg_rot), ' neglected'] );
% disp( ['Second crit: lambda2_ils_tolerance_lath > ' num2str(lambda2_ips_tolerance_lath), ' --> ' num2str(neg_far_from_IPS), ' neglected'] );
% disp( ['Third crit: theta_CP > ', num2str(theta_CP_max), ' degree --> ', num2str(neg_theta_CP) ,' neglected'] );

end