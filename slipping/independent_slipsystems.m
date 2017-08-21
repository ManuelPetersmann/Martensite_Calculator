function [ n_red, d_red ] = independent_slipsystems( plane_fams, dir_fams)
% call:  independent_slipsystems( plane_fams, dir_fams )
% this function takes a Nx3 array of plane family miller index vectors "plane_fams"
% and a Nx3 array of direction family miller index vectors "dir_fams" and 
% returns all combinations of possible slip systems

%% assemble array of full family from family indizes
for i = 1 : size(plane_fams, 1)
    perms = all_from_family_perms( plane_fams(i,:) );
    if i==1
        nn = cat(1, perms );
    else
        nn = cat(1, nn, perms );
    end
end
%
for j = 1 : size(dir_fams, 1)
    perms = all_from_family_perms( dir_fams(j,:) ); 
    if j==1
        dd = cat(1, perms ); 
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
        if abs( dot(m_e, l_e) ) > 1e-5
            % criterion moved from function "reduce_ambiguous_slipsystems" to here
            % 31.7.2017
            continue
        end
        nnn(nr,:) = nn(i,:); 
        ddd(nr,:) = dd(j,:);
        nr = nr +1;
    end
end
%length( nnn )
%length( ddd )

%% reduction to unequivalent systems and those with slip-plane vector
%  to slip direction
[n_red, d_red] = reduce_ambiguous_slipsystems( nnn, ddd );
%
display('Number of individual slip systems is:')
length( n_red )

end


