function [delta_eps, R_dS, new_res, dS] = linesearch_ILS( U, S, S_accummulated, u, ...
                                          delta_eps, initial_res, vec_residual, cpps_gamma )
%
%% Parameters here are important - is there a physical interpretation?
% now i have a trust-region like iteration due to min_delta_res
%min_delta_res = 1.e-6;
dd_eps = delta_eps / 10.;
old_res = initial_res;
%
%disp('left linesearch_while');

% First step
dS = eye(3) + (delta_eps + dd_eps) * S;
[new_res, R_dS] = perp_ILS(S_accummulated*U, dS, u);
%
LT    = R_dS * U;
[ theta_cp_old, ~ ] = min_misorientation( cpps_gamma, LT, true );
%
delta_res = initial_res - new_res;
%delta_theta_cp = 1.;
if ( old_res > new_res )
    
    while   ( old_res > vec_residual )  && ...        % optimization goal is reached
            ( old_res > new_res )
%            ( delta_res / delta_eps > 0.001 ) %&& ...  % improvement is not good enough
    %        ( delta_theta_cp > 0 )                    % MULTI OBJECTIVE OPTIMIZATION 
    %        ( delta_res > vec_residual )              % change of solution is not good enough
        
        delta_eps = delta_eps + dd_eps;
        old_res   = new_res;
        dS = ( eye(3) + (delta_eps + dd_eps) * S );
        [new_res, R_dS] = perp_ILS(U*S_accummulated, dS, u);
        % PET 21.11.17 - multi objective optimization (MOO) for theta_CP --> min
        LT    = R_dS * U;
        [ theta_cp_new, ~ ] = min_misorientation( cpps_gamma, LT, true );
        if(old_res < new_res)
            dd_eps = 0.5 * dd_eps;
        end
        delta_res = old_res - new_res;
        %delta_theta_cp = theta_cp_old - theta_cp_new;
        new_res = new_res * (theta_cp_new / theta_cp_old);
        theta_cp_old = theta_cp_new;
    end
    
end



