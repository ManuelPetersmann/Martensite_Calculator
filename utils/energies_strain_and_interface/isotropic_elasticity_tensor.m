function C = isotropic_elasticity_tensor(K,mu)
%
I = eye(3);
%
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C(i,j,k,l) = C(i,j,k,l) + K * I(i,j) * I(k,l) + mu * (I(i,k)*I(j,l) + I(i,l)*I(j,k) - 2/3* I(i,j) * I(k,l);
            end
        end
    end
end

end

