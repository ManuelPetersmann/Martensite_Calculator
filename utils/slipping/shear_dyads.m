function [ds, ns, S] = shear_dyads(martensite, austenite, miller_dyads)
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

end

