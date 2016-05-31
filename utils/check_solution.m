function [ noSolution, lambda2_smaller1] = check_solution( lambda_1, lambda_2, lambda_3, epsilon )
% This function takes the eigenvalues of a matrix
% and checks if it provides an invariant plane 
% --> eigenvalue lambda_2 ~= 1.0 && lambda_1 < 1.0 && lambda_3 > 1.0

noSolution = true;

% check if eigenvalue lambda_2 is near enough to 1.0 
% --> deviation less than epsilon
if ( abs(lambda_2 - 1.) > epsilon )
  % further check: if deviation greater than epsilon, is lambda_2 greater or
  % smaller than 1.0
    if (lambda_2 - 1.) > 0.
        lambda2_smaller1 = false;
    else
        lambda2_smaller1 = true;
    end
% else - lambda2 solution within precision  - check if the other eigenvalues
% straddle lambda2 = 1, i.e. lambda1 > 1. , lambda3 > 1. 
else if (lambda_1 < 1.)  && (lambda_3 > 1.)
        noSolution = false;
    end
end

end

