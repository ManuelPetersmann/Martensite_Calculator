function funcs = matrix_functions()
funcs.cofactor =                @cofactor;
funcs.polardecomposition =      @polardecomposition;
funcs.is_positive_definite =    @is_positive_definite;
funcs.is_symmetric =            @is_positive_definite;
end

function cofactor( A )
% Calculates the cofactor of matrix A
%TODO
end
%-------------------------------------------------------------------------

function[U,R] polardecomposition( F )
% calculates the polar decomposition of a matrix into its rotational part R
% and its symmetric positive definite part U (pure-distrotion, Bain-strain)
U = sqrt(F'*F) %F'=komplex conjugate transponse!
R = Fi*inv(U0)
end
%-------------------------------------------------------------------------   
% TODO
% function bool =  is_positive_definite( A )
% end
%-------------------------------------------------------------------------
function lol = is_symmetric( A )
if equal_matrices(A,A',1.e-9)
    bool = true;
else
    bool = false;
end
end
%-------------------------------------------------------------------------
