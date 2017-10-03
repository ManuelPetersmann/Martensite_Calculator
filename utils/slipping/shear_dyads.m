function [outputArg1,outputArg2] = shear_dyads(considered_plasticity,ds_product,ns_product,)
%SHEAR_DYADS - TODO - could be integrated like this since more functions
%use the same code below...

martensite.slip_planes, martensite.slip_directions

%% transform product phase slip systems to parent phase and combine all in one array
if considered_plasticity == 2 % only parent phase slip systems
    ds = ds_product;
    ns = ns_product;
end
if considered_plasticity == 1
    for is = 1:size(ds_product,1)
        % transform product phase slip systems to parent phase ones
        ds(is,:) = cp * ds_product(is,:)';
        ns(is,:) = inverse(cp)' * ns_product(is,:)';
    end
end
if considered_plasticity == 3  % if both parent and product phase systems are given
    ds = cat(1,ds,ds_parent);
    ns = cat(1,ns,ns_parent);
    % for outputting found slip systems in miller indizes
    ds_product = cat(1,ds_product,ds_parent);
    ns_product = cat(1,ns_product,ds_parent);
end
%ds
%ns
% for comparison with non-normed vectors...
% for is1 = 1:size(ds,1) % loop over vectors to construct all slip systems
%     S(:,:,is1)  = I + (ds(is1,:)' * ns(is1,:) );
% end
% S(:,:,1)
% S(:,:,80)

% norm vectors and assemble slip systems
for jj = 1:size(ds,1)
    ds(jj,:) = ds(jj,:) / norm(ds(jj,:));
    ns(jj,:) = ns(jj,:) / norm(ns(jj,:));
    S(:,:,jj)  = (ds(jj,:)' * ns(jj,:) );
end

end

