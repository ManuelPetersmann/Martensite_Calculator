 % Note that orientation relationships (ORs) 

function funcs = lattice_correspondance
funcs.correspondance_matrix_components = @correspondance_matrix_components;
funcs.basis_relation =                   @basis_relation;
end

function cp = correspondance_matrix_components(lattice1, lattice2, m1, m2, use_planes)
% this function takes two lattice objects and two matrices of 3 or 2 (for
% orthogonal lattices) direction [...], or plane (...) components (as colums) 
% of each lattice specifying the orientation relationship between the lattices
% i.e. - e.g. m1(:,i) || m22(:,i).  
% The fifth argument is a boolean that specifies wheter planes should be used
% for the OR instead of directions, which are specified by default.
% From this information then the correspondance matrix transforming the first
% set of components to the second is calculated (The reverse transformation
% matrix is given by the inverse). See e.g. Bhadeshia WEGC - p.17
% Principally the Equation: m2 = cp * m1   is solved for cp
if nargin < 5
    use_planes = false;
end
if size(m1,2) == size(m2,2)
    if size(m1,2) == 2
        basfuncs = basis_functions();
        if basfuncs.isorthogonal(lattice1.C) && basfuncs.isorthogonal(lattice2.C)
            m1(:,3) = cross(m1(:,1), m1(:,2));
            m2(:,3) = cross(m2(:,1), m2(:,2));
        else
            error('for non orthogonal lattices three vector pairs are required to calculate the lattice correspondance!')
        end
    end
    % assemble real vectors from components/coordinates
    vec1 = zeros(3);
    vec2 = zeros(3);
    for i = 1:size(m1,2)
        vec1(:,i) = lattice1.vec_from_coords( m1(:,i) );
        vec2(:,i) = lattice2.vec_from_coords( m2(:,i) );
    end
    k = ( norm(vec1(:,1)) / norm(vec2(:,1)) );
    g = ( norm(vec1(:,2)) / norm(vec2(:,2)) );
    m = ( norm(vec1(:,3)) / norm(vec2(:,3)) );
    new = cat(2, k*m2(:,1), g*m2(:,2), m*m2(:,3) );
    %
    % matfuncs = matrix_functions();
    cp = new * inv(m1); % cp Transforms the components/coordinates
    % cp2 = vec2 * matfuncs.inverse(vec1)  This is not the same! only true for basis vectors! 
    % m2 = cp*m1  - is equal to Eq 3.31 Ulrich Müller p.31 (cp = P^-1)
    % Therefore according to 3.29 there cp^-1 transforms the bases into each other
    % Alternatively, instead of inverting cp by calculation, cp^?1 can be deduced by
    % derivation of the matrix for the reverse transformation in the same manner
    % as is done in the second function in this file
    %
    if use_planes
        cp = inverse( cp ); % Transforms components/coordinates
    end
else
    error('the matrices specifying the orientation relationship must be the same size')
end
end
%-------------------------------------------------------------------------
function P = basis_relation( basis1, basis2 )
mf = matrix_functions;
P = basis1 * inv( basis2 );
% Difference between here and Ulrich Müller:
% He writes: basis1 = basis2 * P -> the Basis Vectors are assembled as row
% vectors in the matrices. Here it is written: basis1=P*basis2 -> the Basis
% Vectors are assembled as colum vectors in the matrices. 
% Note (AB)' = B'A'
% Also written as [B1 |> B2] see Cayron 2006 - Groupoid of orientational
% variants
end


% Note: A cross product between two vectors given in the real basis, give a
% vector whose components are given in the reciprocal basis and vice versa - see
% Bhadeshia- Worked examples p.24

%-------------------------------------------------------------------

% based = transform_directions(vectorlist ,vectors, direction )
% transform a list of directions from initial phase
% to final phase based on lattice correspondence

%[ transformed_directions ]
% lattice_correspondence =
% list_of_directions =
% list_of_planes =
% 
% if dircetion == 'toInitial'
% end
%     
% if direction == 'toParent'
%    end
%    
% if vectors == 'planes'
%    end
%    
% if vectors == 'directions'
% end
%   
%-------------------------------------------------------------------


    