function db = calc_dual_basis( basis )
metric = calc_metric( basis );
%mf = matrix_functions;
%inv_metric = mf.inverse(metric);
for i = 1:size( obj.C, 2)
    db(:,i) = inv(metric) * basis(:,i); % Z^i = Z^ij Z_j, does not matter which index
    %is contracted since the metric is symmetric
end
end

