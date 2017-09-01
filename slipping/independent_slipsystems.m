function [ n_red, d_red ] = independent_slipsystems( plane_fams, dir_fams, count_directions_extra)
% call:  independent_slipsystems( plane_fams, dir_fams )
% this function takes Nx3 arrays of plane family miller index vectors "plane_fams"
% and a Nx3 array of direction family miller index vectors "dir_fams" and 
% returns all combinations of possible slip systems,
% where a slip system is defined by a crystallographic plane normal and a
% perpendicular slip direction
% the third (optional) parameter

if nargin < 3
    count_directions_extra = false;
end

%% assemble array of full family from family indizes
for i = 1 : size(plane_fams, 1)
    perms = all_from_family_perms( plane_fams(i,:) );
    if i==1
        nn = perms;
    else
        nn = cat(1, nn, perms );
    end
end
%
for j = 1 : size(dir_fams, 1)
    perms = all_from_family_perms( dir_fams(j,:) ); 
    if j==1
        dd = perms; 
    else
        dd = cat(1, dd, perms );
    end
end

%% assemble all slip systems, i.e. combinations of plane - direction pairings,
% even those not normal to each other
nr = 1;
nnn = [0. 0. 0.];
ddd = [0. 0. 0.];
for i = 1:size(nn,1)
    for j = 1:size(dd,1)
        if abs( dot(nn(i,:),dd(j,:)) ) < 1e-5  % select out non-orthogonal pairings of planes and directions
            nnn(nr,:) = nn(i,:);
            ddd(nr,:) = dd(j,:);
            nr = nr +1;
        end
    end
end
%length( nnn )
%length( ddd )

%% reduction to unequivalent systems and those with slip-plane vector
%  to slip direction
%n_red = nnn;
%d_red = ddd;
[n_red, d_red] = reduce_ambiguous_slipsystems( nnn, ddd, count_directions_extra );
%
display('Number of individual slip systems is:')
length( n_red )

end


