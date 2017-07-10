function[U,R] = polardecomposition( F )
% calculates the polar decomposition of a matrix into its rotational part R
% and its symmetric positive definite part U (pure-distrotion, Bain-strain)
U = sqrtm(F'*F); %F'=komplex conjugate transponse!
R = F*inv(U);
end

