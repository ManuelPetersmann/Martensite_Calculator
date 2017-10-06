function metric = calc_metric( basis )
% call metric( B ) - where B is a matrix containing the basis vectors 
% of a lattice as colums.
metric = zeros(3);
for i=1:3
    for j=1:3
        metric(i,j) = dot( basis(:,i), basis(:,j));
    end
end
end

