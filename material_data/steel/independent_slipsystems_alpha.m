function [ n_red, d_red ] = independent_slipsystems_alpha()
% given the laue symmetry group determines the number of independent slip
% systems of the bcc - martensite system


%% slip systems b.c.c  or  f.c.c
% since the shear is a substantial part of the transformation only 
% shear systems which are favorable in the b.c.c. lattices are considered. 
% the plane and direction families are {110}_alpha, {112}_alpha,
% <111>_alpha, <110>_alpha

%n = cat( 1, all_from_family( [1. 1. 0.],SG ), all_from_family( [1. 1. 2.],SG ) );  
%d = cat( 1, all_from_family( [1. 1. 0.],SG ), all_from_family( [1. 1. 1.],SG ) );

%% or initialisation via permutations of expected indices of miller index families

% get all combinations of indizes
nn = cat(1, perms( [1. 1. 0.] ), perms( [-1. 0. 1.] ) ); % = {1 1 0}perms
nn = cat(1, nn, perms( [1. 1. 2.] ),  perms( [-1. -1. -2.]), perms( [-1. 1. 2.]), perms([-1. 1. -2.]), perms([1. 1. -2.])  ); % = {1 1 2} perms
% add more to get more than 42 systems
nn = cat(1, nn, perms( [-1. -1. 0.] ) );
nn = cat(1, nn, perms( [1. -1. -2.] ), perms( [1. -1. 2.] ), perms( [-1. -1. 2.]) );

dd = cat(1, [1. 1. 1.], [-1. -1. -1.], perms( [1. 1. -1.] ) ); % = <1 1 1>perms
dd = cat(1, dd, perms( [1. 1. 0.] ), perms( [-1. 0. 1.] ) ); % = <1 1 0>perms
% add more to get more than 42 systems
dd = cat(1, dd, perms( [-1. -1. 1.] ), perms( [-1. 1. 1.] ) );
dd = cat(1, dd, perms( [-1. -1. 0.] ), perms( [1. 0. -1.] ) );

nn = reduce_vecs(nn);
dd = reduce_vecs(dd);

% get all combinations of plane - direction pairings, even those not normal
% to each other
nr = 1;
nnn = [0. 0. 0.];
ddd = [0. 0. 0.];
for i = 1:size(nn,1)
    for j = 1:size(dd,1)
        nnn(nr,:) = nn(i,:); 
        ddd(nr,:) = dd(j,:);
        nr = nr +1;
    end
end

%% reduction to unequivalent systems 
[n_red, d_red] = reduce_ambiguous_slipsystems( nnn, ddd );
length( n_red )
length( d_red )

%cat(2, n_red, d_red)


end


