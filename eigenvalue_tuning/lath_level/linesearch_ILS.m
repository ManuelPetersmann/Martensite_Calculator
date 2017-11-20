function [delta_eps, R_dS, new_res, dS] = linesearch_ILS( U, S, S_accummulated, u, ...
                                                          delta_eps, initial_res, vec_residual )
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
[new_res, R_dS] = perp_ILS(U, dS, S_accummulated, u);
delta_res = old_res - new_res;
if ( old_res > new_res )
    
    while   ( old_res > vec_residual )  && ...        % optimization goal is reached
            ( delta_res / delta_eps > 0.001 ) %&& ...  % improvement is not good enough
    %        ( delta_res > vec_residual )              % change of solution is to
        
        delta_eps = delta_eps + dd_eps;
        old_res   = new_res;
        
        dS = ( eye(3) + (delta_eps + dd_eps) * S );
        [new_res, R_dS] = perp_ILS(U, dS, S_accummulated, u);
        if(old_res < new_res)
            dd_eps = 0.5 * dd_eps;
        end
        delta_res = old_res - new_res;
        
    end
    
end
