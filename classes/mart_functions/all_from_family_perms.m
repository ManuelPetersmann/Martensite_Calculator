function [ all ] = all_from_family_perms( vec, sort_out_negatives )
% this function takes a vector of three miller indizes representing 
% a plane or direction family and returns the full family
% "negatives" is a boolean specifying if sign-ambiguous entries [multiplied
% by (-1)] should be sorted out (true) e.g. because the vector [1 2 3] equals [-1 -2 -3]

if nargin < 2
    sort_out_negatives = false; % --- % edited 29.08.17 Manuel
end

x = vec(1);
y = vec(2);
z = vec(3);
ambiguous = cat(1,            perms( [x y z]   ), perms( [-x y z] ),  perms( [x -y z] ),  perms( [x y -z] )   ); 
ambiguous = cat(1, ambiguous, perms( [-x -y z] ), perms( [x -y -z] ), perms( [-x y -z] ), perms( [-x -y -z] ) );

all = reduce_vecs( ambiguous, sort_out_negatives );  


end

