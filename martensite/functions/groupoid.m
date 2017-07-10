function [ gpt, operators ] = groupoid( T_aust_mart_i )       

% TODO TEST and add "gpt" output matrix

% The set of variants (o_ij) and operators form the groupoid
% of orientational variants (see Cayron 2006 Acta Crystall)
% o_ij is a set of isometry matrices that contains at least one matrix with positive
% determinant, i.e. a rotational matrix
%
% This function takes the point group "pg" (3x3xn) Matrix and the Bain
% strain "B" and calculates the 1) GROUPOID TABLE
%                 o1 o2 ..... o_n  allg U1->Ui
%                 variants...
% o1 , variant   [o_ij] = (o_m)^-1 * o_n   m...colum operator, n...row operator
% o2 , variant
% ... allg Ui -> U1

% The product of operator [o_ij] is generally multivalued e.g. two mappings
% map a1 to a2 and three mappings map a1 to a3 --> o_13 = 2x3 matrix (with
% a sum of 6 values)

% and ii) the OPERATOR LIST "operators" with entries
% [ operator nr, angle,   0   1 4 ...                                    ]
% [ axis_comp1,    ac2, ac3   1 2 ...pairs of operators [1 1] [4 2] etc. ] 
% ----------------------------------------------

mf = matrix_functions;
operators = zeros(2,3);

for i=1:size(T_aust_mart_i,3)
    for j=i+1:size(T_aust_mart_i,3)
        A =  mf.inverse( T_aust_mart_i(:,:,i) ) * T_aust_mart_i(:,:,j);
        % solve characteristic equation to find the rotation axis
        % (corresponding eigenvalue is 1)
        [V,D] = eig(A);
        % D... Diagonal matrix of eigenvalues
        % V... colum vectors are correspondig right eigenvectors i.e. A*V = V*D
        if ( abs(det(A) - 1) < 1.e-4 )
            for k=1:size(V,2)
                if ( abs(V(k,k) - 1) < 1.e-4 )
                    axis = D(:,k);
                    angle = acos((trace(A)-1.)/2.);
                    break
                end
            end
            [is_new, nr] = new_operator(axis, angle);
            if is_new
                operators(:,:,size(operators,3)+1) = [size(operator,3)+1 , angle, 0; 
                                                    axis'];
            else 
                operators(:, size(operators,2)+1 ,nr) = [i,j];
            end            
        end
    end
end

% nested function
    function [bool, nr] = new_operator(ax, ang)
        bool = true;
        nr = -1;
        for m = 1:size(operators,3)
            if ( abs(operators(2,1:3,m) - ax) < 1.e-4) && ( abs(operators(1,2,m) - ang) < 1.e-4 )
                bool = false;
                nr = m;
                break
            end
        end
    end
    
end % end main function


