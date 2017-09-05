function [ isSolution, lambda2_smaller1] = check_IPS_solution( lambda_1, lambda_2, lambda_3, epsilon )
% This function takes the eigenvalues of a matrix in inreasing order
% and checks if the matrix depicts an invariant plane mapping 
% --> eigenvalue lambda_2 == 1.0 && lambda_1 < 1.0 && lambda_3 > 1.0
% as a forth optional argument the tolerable deviation lambda2 -1. = epsilon 
% can be specified. The default is 1.e-9 

if nargin < 4
    epsilon = 1.e-12;
end

isSolution = false;

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
elseif ( (lambda_1 < 1.)  && (lambda_3 > 1.) )
        isSolution = true;
        lambda2_smaller1 = 999999.; % assign any value so no error occurs
end

end

