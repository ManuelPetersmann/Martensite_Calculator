function [U,R] = polardecomposition( F )
% call: [U,R] = polardecomposition( F )
% calculates the polar decomposition of a matrix into its rotational part R
% and its symmetric positive definite part U (pure-distrotion, Bain-strain)
% matlab uses the Higham numerical method for 3x3 matrices !!! see documentation 
U = sqrtm(F'*F); %F'=komplex conjugate transponse!
R = F*inverse(U);
end

