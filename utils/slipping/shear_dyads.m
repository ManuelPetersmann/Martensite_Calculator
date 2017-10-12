function [ds, ns, S, slip_combinations] = shear_dyads(martensite, austenite, miller_dyads)
% SHEAR_DYADS - this function takes information on the slip systems
% and creates the necessary matrices used in the middle eigenvalue modification function
% the first two inputs are the martensite and austenite objects, the third
% input determines wheter the vectors forming the shear dyads are normed or
% if the miller indizes are used directly

if nargin < 3
    miller_dyads = false;
end

%% transform product phase slip systems to parent phase and combine all in one array
if (martensite.considered_plasticity == 1 || martensite.considered_plasticity == 3) % only product/martensite phase slip systems
    for is = 1:size(martensite.slip_directions,1)
        % transform product phase slip systems to parent phase ones
        ds(is,:) = martensite.cp * martensite.slip_directions(is,:)';
        ns(is,:) = inverse(martensite.cp)' * martensite.slip_planes(is,:)';
    end
end
%
if martensite.considered_plasticity == 2 % only parent/austenite phase slip systems
    ds = austenite.slip_directions;
    ns = austenite.slip_planes;
end
%
if martensite.considered_plasticity == 3  % if both parent and product phase systems are given
    ds = cat(1,ds,austenite.slip_directions);
    ns = cat(1,ns,austenite.slip_planes);
end

% norm vectors and assemble slip systems
for jj = 1:size(ds,1)
    if ~miller_dyads
        S(:,:,jj) = ( ds(jj,:) / norm(ds(jj,:)) )' * ( ns(jj,:) / norm(ns(jj,:)) );
    else
        S(:,:,jj)  = (ds(jj,:)' * ns(jj,:) );
    end
end

%% give the martensite slip systems in the martensite basis to use miller indizes for notion
% make 1x4 vectors with the 4-th entry specifying wheter the system is from
% martensite 99 or austenite 88 - also useful for writing results
if (martensite.considered_plasticity == 1 || martensite.considered_plasticity == 3)
    % append identification NR (99) to austenite slip systems --> [s,99], [n,99] e.g. [1 1 1 99]
    phase_identifier = 99*ones(size(martensite.slip_directions,1),1);
    ds = cat(2,martensite.slip_directions,phase_identifier);
    ns = cat(2,martensite.slip_planes,    phase_identifier);
end
%
if (martensite.considered_plasticity == 2 || martensite.considered_plasticity == 3)
    % append identification NR (88) to austenite slip systems --> [s,88], [n,88] e.g. [1 1 1 88]
    phase_identifier = 88*ones(size(austenite.slip_directions,1),1);
    ds_aust = cat(2,austenite.slip_directions,phase_identifier);
    ns_aust = cat(2,austenite.slip_planes,    phase_identifier);
    if martensite.considered_plasticity == 2
        ds = ds_aust;
        ns = ns_aust;
    end
end
%
if martensite.considered_plasticity == 3
    ds = cat(1,ds,ds_aust);
    ns = cat(1,ns,ns_aust);
end

slip_combinations = nchoosek(size(ds,1),2);
disp( ['Number of possible pairings is = ', num2str( slip_combinations )])
disp('nr of solutions cannot be greater than 2-times this value.')

end

