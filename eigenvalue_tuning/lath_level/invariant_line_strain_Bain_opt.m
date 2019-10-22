function solutions = invariant_line_strain_Bain_opt(martensite, austenite)
% calculates invariant line strain (rototated, unstreched) - per default
% the invariant line is assummed to be the close packed direction in
% austenite!

if isempty( martensite.invariant_lines )
    us = austenite.CP_dirs
    % VECTORS MUST BE NORMED!
end

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
	theta_cp_old = 0.;
    %
    for is1 = 1:length(S)
        for is2 = is1+1:length(S)
            
            % res_initial = norm( u - u2_unrot ); % If an invariant line is specified a -
            % res_old = res_initial;
            
            %% all this is calculated only to check if an improvement is possible to begin with - 
            delta_eps = delta_eps_initial;
            eps1 = eps_initial;
            eps2 = eps_initial;
            S_accummulated = eye(3);
            %
            R_dS1 = eye(3);
            R_dS2 = eye(3);
            shear_increments = [];

            while ( res_old > vec_residual )  ...
                    && (eps1 < eps_max)        && (eps2 < eps_max)
                
                % if the solution for 'eps' is very high or low respectively, do not consider it
                if ((eps1 < eps_initial) || (eps2 < eps_initial))
                    error('this should not happen - fix code...')
                end
                
                % dS1 = eye(3) + delta_eps*S(:,:,is1);
                % dS2 = eye(3) + delta_eps*S(:,:,is2);
                dS1 = R_dS1 * ( eye(3) + delta_eps*S(:,:,is1) ) * R_dS1';
                dS2 = R_dS2 * ( eye(3) + delta_eps*S(:,:,is2) ) * R_dS2';
                
                [res_B1, res_vec1, RS1, theta_cp_new1] = bain_opt_shear( martensite.U, dS1, S_accummulated, u )
                [res_B2, res_vec2, RS2, theta_cp_new2] = bain_opt_shear( martensite.U, dS2, S_accummulated, u )
                
                LT1    = R_dS1 * martensite.U;
                [ theta_cp_new1, ~ ] = min_misorientation( austenite.CPPs, LT1, true );
                
                new_optfunc1 = new_res_dS1 * (theta_cp_new1 / theta_cp_old);
                

                %if(new_res_dS1 < new_res_dS2)
                if(new_optfunc1 < new_optfunc2)
                    % here generally the minimum should be taken if the solution has not been passed already
                    % from one search direction, otherwise the shear amplitude is adopted first and only then it is
                    % checked again wheter a new minimum does not overshoot the solution from one direction
                    %if ( res_old < new_res_dS1 )
                    if ( old_optfunc < new_optfunc1 )
                        % i.e. the solution is passed or delta_lambda2_to_1 does not decrease, which should not happen - cut back
                        delta_eps = stepwidth_change * delta_eps;
                    else
                        old_optfunc = new_optfunc1;
                        theta_cp_old = theta_cp_new1;
                        res_old = new_res_dS1;
                        eps1 = eps1 + delta_eps;
                        S_accummulated = S_accummulated * dS1;
                        R_mod = R_dS1;
                        shear_increments( : , size(shear_increments,2)+1 ) = [delta_eps; 0];
                    end
                else
                    %if ( res_old < new_res_dS2 )
                    if ( old_optfunc < new_optfunc2 )
                        delta_eps = stepwidth_change * delta_eps;
                    else
                        old_optfunc = new_optfunc2;
                        theta_cp_old = theta_cp_new2;
                        res_old = new_res_dS2;
                        eps2 = eps2 + delta_eps;
                        S_accummulated = S_accummulated * dS2;
                        R_mod = R_dS2;
                        shear_increments( : , size(shear_increments,2)+1 ) = [0; delta_eps];
                    end
                end
                
                if (delta_eps < delta_eps_tolerance) % 0.5^16 = 1.e-5||
                    break
                end
                
                
            end % end while
            

            if res_old < 0.005
                res_old
                theta_cp_old
            end
            
            if ( res_old < vec_residual ) && (theta_cp_old < 5.) 
                
                F_tot = R_mod * martensite.U *  S_accummulated;
                LT    = R_mod * martensite.U;     
                
                isol = isol + 1;
                %if mod(isol,10)==0
                    isol
                    %pause(1);
                %end
                
                eps_s = [eps1; eps2];
                d = [ds(is1,:); ds(is2,:)];
                n = [ns(is1,:); ns(is2,:)];
                
                %solutions.array( isol-1 ) =  ILS_solution(u, F_tot, LT); % F_tot = ST...shape transformation
                %solutions.array( isol-1 ).slip = Slip_systems( eps_s, d, n );
                solutions.array( isol )      =  ILS_solution(u, F_tot, LT, R_Bain);
                solutions.array( isol ).shear_increments = shear_increments;
                solutions.array( isol ).id   = isol;
                solutions.array( isol ).slip = Slip_systems( eps_s, d, n );    
            %else
            %    neg_no_convergence_to_ILS = neg_no_convergence_to_ILS +1;
            end
                
        end
    end
end
comb = nchoosek(length(S),2)*size(us,1);
 disp( [num2str(comb),' combinations {u_i,S1_j,S2_k} tested, ', ...
     num2str(isol),' solutions found, ' num2str(comb-isol), ' did not converge to ILS'] );
% disp( ['First crit: rotation target function of lath > ', num2str(max_rot_angle_inclusion),' degree --> ', num2str(neg_rot), ' neglected'] );
% disp( ['Second crit: lambda2_ils_tolerance_lath > ' num2str(lambda2_ips_tolerance_lath), ' --> ' num2str(neg_far_from_IPS), ' neglected'] );
% disp( ['Third crit: theta_CP > ', num2str(theta_CP_max), ' degree --> ', num2str(neg_theta_CP) ,' neglected'] );


% FSF_bhadeshia = [1.125317  -0.039880   0.106601;
%                  0.023973   1.129478   0.084736;
%                 -0.154089  -0.115519   0.791698]



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

