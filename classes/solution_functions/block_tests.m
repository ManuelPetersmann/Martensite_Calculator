ffunction [block_sols, block_solutions] = block_tests(lath_solutions, block_solutions, U, lambda2_tol, cof_tol, det_tol) % outarg - block_solutions
% call: block_tests(lath_solutions, block_solutions,tol)
% function to investigate solution space of blocks
%
% lath_solutions ... array of lath solutions for building blocks
% block_solutions ... object of Solution_array_composite a.o. with property

% F_composite = x*F_lath1 + (1-x)*F_lath2
% with x being the volume fraction of solution 1.
% Note, that for such an average the IPS condition is approximately
% maintained, which has been checked to be true for various combinations.

% Note: THE AVERAGE OF TWO IPS ALWAYS IS AN IPS (except for the small 
% error that is made because the determinant is not invariant to addition)
% In general the linear mixture rule is only valid in the case of
% n->infinity (minor relations) see Bhattacharya - Microstructures of
% martensties - p.131.dö

if nargin < 4
    lambda2_tol_block_aust = 1.e-3 % doesnrt matter if 0.001 or 0.0001 !!! important! some more solutions with 0.003
    lambda2_tol_laths = 1.e-2
    cof_tol = 1.e-2
    det_tol = 1.e-2
    block_hp_cp_aust_tol = 5.; % degree - even if i just set this only to 10 most solutions fall out  
    %
    rot_angle_block = 3.
end

    function lambda2_mix = mix_y2( x, F1, F2)
        Fc = linmix2(x, F1, F2);
        [~,lambda2_mix] = sorted_eig_vals_and_vecs(Fc'*Fc);
    end

detU = det(U);
I = eye(3); % = austenite
count = 0;
block_sols = 0;


neg_minors = 0;
%
neg_rot_angle = 0;
%
neg_lamda2_block_aust = 0;
neg_lamda2_laths = 0;
neg_hp = 0;

%neg_mix_res_hc = 0;
%neg_mix_res_hh = 0;

%xi = linspace(0,1,100);
%xi = linspace(0.1,0.85,6);
xi = linspace(0,1,10);


% loop over slip system combinations
for is1 = 1: (size(lath_solutions.array,2)-1)
    for is2 = (is1+1): size(lath_solutions.array,2)
        sol1 = lath_solutions.array(is1);
        sol2 = lath_solutions.array(is2);    
        F1 = sol1.ST;
        F2 = sol2.ST;
        
        count = count + 1;
        
        x = 0.5;
        Fc = linmix2(x,F1,F2);
        
                cofFc = cofactor( Fc );
        cof_F_sum = x * cofactor(F1)  +  (1.-x) * cofactor(F2);
        sum(sum(abs(cofFc - cof_F_sum)))
     
        
        %% third MINORS RULE
        det_Fc = det( Fc ); % plotting showed that if the determinant changes then the maximum deviation is at xi=0.5
        if abs(detU - det_Fc) > det_tol
            neg_minors = neg_minors + 1;
            continue
        end
%         
        %% second MINORS RULE
        cofFc = cofactor( Fc );
        cof_F_sum = x * cofactor(F1)  +  (1.-x) * cofactor(F2);
        if sum(sum(abs(cofFc - cof_F_sum))) > cof_tol % frob_distance(cofFc , cof_F_sum)
            neg_minors = neg_minors +1;
            continue
        end
        
        
        [~,R] = polardecomposition( Fc );
        vec4 = vrrotmat2vec( R );
        % convert ange to degree
        angle = rad2deg( vec4(4) );
        if angle > rot_angle_block
            neg_rot_angle = neg_rot_angle +1;
            continue
        end
%        angle

                
%         %% RANK one between block-aust -check deviation of lambda2
        if (mix_y2(x,F1,F2) - 1)  >  lambda2_tol_block_aust
            neg_lamda2_block_aust = neg_lamda2_block_aust +1;
            % mix_y2(x,F1,F2) - 1
            continue
        end
        
        
        %% RANK one between laths
%         if ~is_rank_one_connected(F1,F2,lambda2_tol_laths)
%             neg_lamda2_laths = neg_lamda2_laths + 1;
%             continue
%         end
        
        %% deviation of average block habit plane form 111_aust -- should be sorted out afterwards !!!!
       [y1, y3, d1, d2, h1, h2, Q1, Q2] = rank_one(Fc, I, lambda2_tol_block_aust, false); % last 'false' is that no lambda_2_warning occurs
%        PET 5.11. - corrected Error > instead of < !
        if ( (min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h1) > block_hp_cp_aust_tol) && ... % should here be an && ?
             (min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h2) > block_hp_cp_aust_tol) )
            neg_hp = neg_hp +1;
            continue
        end
%         min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h1)
%         min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h2)
        
%         % Considering that the two F could be equal since the slip
%         deformations are not linearly independent! c.f. non-uniqueness of
%         plastic slip
%         if sum(sum(abs(F1 - F2 ))) < 0.1 %delta_F_min
%         if ( sum(sum(abs( inverse(F1)*(F2) ))) - 3) < delta_F_min
%             neg_diff = neg_diff + 1;
%             continue
%         end

        cofFc = cofactor( Fc );
        cof_F_sum = x * cofactor(F1)  +  (1.-x) * cofactor(F2);
        sum(sum(abs(cofFc - cof_F_sum)))

        
        
%         %% deviation of average block habit plane form 111_aust -- should be sorted out afterwards !!!!
%        [y1, y3, d1, d2, h1, h2, Q1, Q2] = rank_one(Fc, I, lambda2_tol, false); % last 'false' is that no lambda_2_warning occurs
% %       min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h1) 
% %       min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h2)
% %        PET 5.11. - corrected Error > instead of < !
%         if ( (min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h1) > block_hp_cp_aust_tol) || ... % should here be an && ?
%              (min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h2) > block_hp_cp_aust_tol) )
%             neg_hp = neg_hp +1;
%             continue
%         end
        
%        lath_solutions.cryst_fams('KS')
%         sol1.id
%         sol2.id        
        
        % do not mix variants not fullfilling predefined criteria
        % e.g. habit plane deviation from {111}_aust or something else
        % up to now there are two criteria ( both angle tolerances )
        % first angle between common line of invariant habit planes and
        % preferred invariant line (e.g. of set of cp-directions)
        % second angle between habit planes       
%         if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
%             % considering longitudinal dimension of lath -> a
%             vec_in_both_planes = cross( sol1.h , sol2.h );
%             % check if habit planes are not parallel
%             if abs(vec_in_both_planes) < 1.e-8 % entry wise for all entries
%                 theta_intersec_cpdir = misorientation_vector_and_plane( lath_solutions.cryst_fams('KS'), sol1.h );
%             else
%                 theta_intersec_cpdir = min_misorientation( lath_solutions.cryst_fams('KS'), vec_in_both_planes );
%             end
%             %
%             if theta_intersec_cpdir  >  block_solutions.mixing_tolerances('theta_intersec_cpdir')
%                 neg_mix_res_hc = neg_mix_res_hc +1;
%                 continue
%             end
%         end
      
%         if isKey(block_solutions.mixing_tolerances,'theta_hps')
%             theta_hps = get_angle( sol1.h , sol2.h );
%             if theta_hps  >  block_solutions.mixing_tolerances('theta_hps')
%                 % considering width of laths -> b
%                 neg_mix_res_hh = neg_mix_res_hh +1;
%                 continue
%             end
%         end   


%         nf11 = norm(sol1.F1(:,1));
%         nf12 = norm(sol1.F1(:,2));
%         nf13 = norm(sol1.F1(:,3));
%         nf21 = norm(sol2.F1(:,1));
%         nf22 = norm(sol2.F1(:,2));
%         nf23 = norm(sol2.F1(:,3));
%         acosd( dot(sol1.F1(:,1) / nf11 , sol2.F1(:,1) / nf21 ) )
%         acosd( dot(sol1.F1(:,2) / nf12 , sol2.F1(:,2) / nf22 ) )
%         acosd( dot(sol1.F1(:,3) / nf13 , sol2.F1(:,3) / nf23 ) ) 
%         nf11 - nf21
%         nf12 - nf22
%         nf13 - nf23
% 
% %        sum(sum(abs(F1 - F2 )))
% %        sum(sum(abs( inverse(F1)*F2 )))
% %         [~,R] = polardecomposition ( inverse(F1)*F2 );
% %         r = angle_from_rotmatrix( R )
%             F1
%             F2
%             [sol1.shear_direction, sol1.slip_normal_plane_vec]
%             [sol2.shear_direction, sol2.slip_normal_plane_vec]
%             sol1.eps_s
%             sol2.eps_s
%             sol1.stepwidth
%             sol2.stepwidth


        
%          for i=1:length(xi)
%              x(i) = xi(i);
% % %             y2(i) = mix_y2( x, F1, F2);
%              Fc = linmix2(x,F1,F2);
% % %             detFc(i) = detU - det(Fc);
% % %             cofFc = cofactor( Fc );
% % %             cof_F_sum = x * cofactor(F1)  +  (1.-x) * cofactor(F2);
% % %             frob_dist(i) = frob_distance(cofFc , cof_F_sum);
% %             try
%                  [y1, y3, d1, d2, h1, h2, Q1, Q2] = rank_one(Fc, I, lambda2_tol, false); % last 'false' is that no lambda_2_warning occurs
%                 hx(i) = h1(1);
%                 hy(i) = h1(2);
%                 hz(i) = h1(3);
% %                 sqrt( h1(1)^2 + h1(2)^2 + h1(3)^2 )
% %             catch
% %                 continue
% %             end
% %            [ theta1(i), closest_from_vecs1(:,i)] = min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h1) %, plane )
% %            [ theta2(i), closest_from_vecs2(:,i)] = min_misorientation( lath_solutions.cryst_fams('cpps_gamma'), h2) %, plane )
%             eps(i) = sqrt(y3) - sqrt(y1);
%         end
%         
%        
%          figure;
%           plot(x,hx,'-o',x,hy,'-o',x,hz,'-o');
          
%         plot(xi,theta1,'-o');
%         hold on
%         plot(xi,theta2,'-o');
       
%         plot(xi,eps,'-o');
%         plot(xi,detFc);
%         legend('delta_det','-o');
%         hold on
%         plot(xi,y2-1.,'-o');
%         plot(xi,frob_dist,'-o')    

%abs(detU - det_Fc)
         block_sols = block_sols + 1;
                   
        block_solutions.array( block_sols ).lath_solution_pair = [sol1, sol2];  % U,tolerance]; %
         
         if mod(block_sols,100)==0
             block_sols
             count
         end         


    end % end of loop 1
end % end of loop 2


count
disp( ['First  crit: ', num2str(neg_minors), ' neglected due to minors relations (cof_tol =',num2str(cof_tol), ', det_tol =',num2str(det_tol),')'] );
neg_rot_angle
disp( ['Second crit: ', num2str(neg_lamda2_block_aust), ' neglected because lamda2 of block deviates more than ', num2str(lambda2_tol_block_aust), ' from 1'] );
disp( ['Third crit: ', num2str(neg_lamda2_laths), ' neglected laths deformations of pairings are not rank one connected with tolerance ', num2str(lambda2_tol_laths) ] ); 
disp( ['Fourth crit: ', num2str(neg_hp), ' neglected because ave HP deviates more than ', num2str(block_hp_cp_aust_tol),' from 111_aust']);
%disp( ['Fourth crit: ', num2str(neg_diff), ' neglected because Fs do not differ in the d1 norm more than ', num2str(delta_F_min) ] ); 

if isKey(block_solutions.mixing_tolerances,'theta_hps')
    disp( [ num2str(neg_mix_res_hh), ' mixings neglected due habit plane angle < ', ...
        num2str(block_solutions.mixing_tolerances('theta_hps') ), ' on mixing laths to blocks'] )
end
if isKey(block_solutions.mixing_tolerances,'theta_intersec_cpdir')
    disp( [ num2str(neg_mix_res_hc), ' mixings neglected due ang(h1 x h2, <110>_aust ) < ',...
        num2str(block_solutions.mixing_tolerances('theta_intersec_cpdir') ) ])
end

disp( ['number of potential solutions found: n_sol = ', num2str(block_sols) ] )




end