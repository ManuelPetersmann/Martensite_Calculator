function [ isSolution, lambda2_smaller1] = check_solution( lambda_1, lambda_2, lambda_3, epsilon )
% This function takes the eigenvalues of a matrix
% and adopts the shear so that the second eigenvalue gets closer to 1.

isSolution = false;

if ( abs(lambda_2 - 1.) > epsilon )
    if (lambda_2 - 1.) > 0.
        lambda2_smaller1 = false;
    else
        lambda2_smaller1 = true;
    end
% else - lambda2 solution within precision  - check if the other eigenvalues
% straddle lambda2 = 1, i.e. lambda1 > 1. , lambda3 > 1. 
elseif ( (lambda_1 < 1.)  && (lambda_3 > 1.) )
        isSolution = true;
        lambda2_smaller1 = 5.; % assign any value so no error occurs
end

end

