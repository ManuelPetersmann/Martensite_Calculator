function [R, alpha] = optimize_rotation_around_invariant_line(inv_vec , BS1S2 , eps)
% call:  optimize_rotation_around_invariant_line(inv_vec,BS1S2)
% the invariant line direction vector (inv_vec) simultaneously is the rotation axis
% R [axis=inv_vec, angle=alpha] * B*S1*S2 x = x (lambda = 1) --> Invariant line strain

if nargin < 3
    eps = 1.e-9;
end
inv_vec = inv_vec / norm( inv_vec );
I = eye(3);

if size(inv_vec) == [1,3]
    inv_vec = inv_vec';
end

%% optimization parameters
res_old = norm( (BS1S2 - I)*inv_vec );
search_direction = 1.;
alpha = deg2rad(50.);
dalpha = deg2rad(1.);

%% determine sign (search direction) of alpha and wheter a rotation can improve the ILS characteristics in the first place 
R_test = axis_angle_to_rotmat( alpha, inv_vec );
new_vec = (R_test*BS1S2 - I)*inv_vec;
res_new = norm(new_vec)
if res_new > res_old
    search_direction = -1;
    dalpha = search_direction * angle_increment;
    % reverse rotation
    alpha = search_direction*alpha;
    R_test = axis_angle_to_rotmat( alpha, inv_vec );
    new_vec = (R_test*BS1S2 - I)*inv_vec;
    res_new = norm(new_vec)
    if res_new > res_old
        error('a rotation about the specified axis cannot improve the overall deformation to be an Invariant Line strain (ILS)')
    else
        delta_res = res_old - res_new
        res_old = res_new;
    end
else
    delta_res = res_old - res_new
    res_old = res_new;
end


%%
counter = 0;
while (delta_res > 1.e-6)  &&  (res_new > eps) && (counter < 1000)  % same eps for tolerance to end iteration and goal! 
    
    R = axis_angle_to_rotmat( alpha + dalpha , inv_vec );
    1. - det(R)
    new_vec = (R*BS1S2 - I)*inv_vec;
    res_new = norm(new_vec)
    
    if (res_new > res_old)% residual grows -> make increment size smaller
        dalpha = dalpha* 0.5;
    else
        % calc improvement and update residual
        delta_res = res_old - res_new
        res_old = res_new;
        %
        % update alpha
        alpha = alpha + dalpha;
        rad2deg(alpha)
    end
    counter = counter + 1
end

end

