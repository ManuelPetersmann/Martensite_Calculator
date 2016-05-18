function metric = metric( basis )
metric = zeros(3);
for i=1:3
    for j=1:3
        metric(i,j) = dot( basis(:,i), basis(:,j));
    end
end
end

